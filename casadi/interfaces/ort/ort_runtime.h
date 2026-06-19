//
//    MIT No Attribution
//
//    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//


/* Runtime support for the black-box ONNX Runtime interface. Shared by the plugin's eval()
 * and by generated C code; both must compile with the ONNX Runtime C API on the include path
 * and link libonnxruntime. Inputs/outputs are double (casadi_real), indices long long.
 * Every numeric ONNX tensor type is supported by converting to/from double. */

#ifndef CASADI_ORT_RUNTIME_H
#define CASADI_ORT_RUNTIME_H

#include <onnxruntime_c_api.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Model metadata; dims are flattened, with in_ndim/out_ndim giving the per-tensor counts.
 * ORT requires every model input to be fed, so n_in/inputs cover ALL model inputs; in_src[i]
 * selects what feeds input i: >=0 a caller arg index, -1 the default (NaN/0), -2 a baked
 * constant taken from in_val (flat over all inputs, numel each). outputs cover only the
 * selected (exposed) outputs. */
struct casadi_onnxruntime_prob {
  long long n_in, n_out;
  const char** input_names;
  const char** output_names;
  const long long* in_src;
  const double* in_val;
  const long long* in_elem_type;
  const long long* out_elem_type;
  const long long* in_ndim;
  const long long* out_ndim;
  const long long* in_dims;
  const long long* out_dims;
  const long long* in_numel;
  const long long* out_numel;
  const unsigned char* model_data;
  long long model_size;
};

/* Per-instance ONNX Runtime state. The session/env/mem handles are created once by
 * casadi_onnxruntime_init (the "prepare" step of the ONNX backend contract). The remaining
 * fields hold the value-independent per-eval scaffolding built once by casadi_onnxruntime_prepare
 * and reused across every casadi_onnxruntime_solve, so repeated evaluation does no allocation:
 *  - buf[i]  : persistent typed input buffer (the OrtValue inv[i] is a non-owning view of it)
 *  - inv[i]  : input tensor, created once over buf[i]; solve() only overwrites buf[i] contents
 *  - outv    : scratch handles owned by Run() and released each solve()
 *  - row     : conversion scratch sized to the largest input/output tensor */
struct casadi_onnxruntime_data {
  OrtSession* session;
  OrtEnv* env;
  OrtMemoryInfo* mem;
  OrtValue** inv;
  OrtValue** outv;
  void** buf;
  double* row;
  int prepared;
};

static const OrtApi* casadi_onnxruntime_api(void) {
  return OrtGetApiBase()->GetApi(ORT_API_VERSION);
}

/* Floating-point ONNX types can hold a NaN sentinel; integer types cannot */
static int casadi_onnxruntime_is_float(long long et) {
  return et == ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT
      || et == ONNX_TENSOR_ELEMENT_DATA_TYPE_DOUBLE
      || et == ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT16
      || et == ONNX_TENSOR_ELEMENT_DATA_TYPE_BFLOAT16;
}

/* IEEE half / bfloat16 <-> double (stored as uint16_t) */
static double casadi_onnxruntime_half2d(uint16_t h) {
  uint32_t sign = ((uint32_t) (h & 0x8000u)) << 16;
  uint32_t e = (uint32_t)((h >> 10) & 0x1f), m = (uint32_t)(h & 0x3ff), bits;
  float f;
  if (e == 0) {
    if (m == 0) { bits = sign; }                       /* signed zero */
    else {                                             /* subnormal -> normalize */
      e = 127 - 15 + 1;
      while ((m & 0x400u) == 0) { m <<= 1; --e; }
      m &= 0x3ffu;
      bits = sign | (e << 23) | (m << 13);
    }
  } else if (e == 31) {
    bits = sign | 0x7f800000u | (m << 13);             /* inf / nan */
  } else {
    bits = sign | ((e + 127 - 15) << 23) | (m << 13);  /* normal */
  }
  memcpy(&f, &bits, 4);
  return (double) f;
}
static uint16_t casadi_onnxruntime_d2half(double d) {
  float f = (float)d; uint32_t b; uint16_t s, e; int32_t ee;
  memcpy(&b, &f, 4);
  s = (uint16_t)((b >> 16) & 0x8000u);
  ee = (int32_t)((b >> 23) & 0xff) - 127 + 15;
  if (ee <= 0) return s;                       /* underflow -> signed zero */
  if (ee >= 31) return (uint16_t)(s | 0x7c00u);/* overflow -> inf */
  e = (uint16_t)ee;
  return (uint16_t)(s | (e << 10) | (uint16_t)((b >> 13) & 0x3ffu));
}
static double casadi_onnxruntime_bf162d(uint16_t v) {
  uint32_t u = ((uint32_t)v) << 16; float f; memcpy(&f, &u, 4); return (double)f;
}
static uint16_t casadi_onnxruntime_d2bf16(double d) {
  float f = (float)d; uint32_t u; memcpy(&u, &f, 4); return (uint16_t)(u >> 16);
}

/* CasADi matrices are column-major, ONNX tensors row-major: rank-2 differs by a transpose */
static void casadi_onnxruntime_to_row(long long ndim, const long long* dims,
                                      const double* col, double* row, long long n) {
  if (ndim == 2) {
    long long r, c, d0 = dims[0], d1 = dims[1];
    for (r = 0; r < d0; ++r) for (c = 0; c < d1; ++c) row[r * d1 + c] = col[r + c * d0];
  } else {
    long long j; for (j = 0; j < n; ++j) row[j] = col[j];
  }
}

static void casadi_onnxruntime_from_row(long long ndim, const long long* dims,
                                        const double* row, double* col, long long n) {
  if (ndim == 2) {
    long long r, c, d0 = dims[0], d1 = dims[1];
    for (r = 0; r < d0; ++r) for (c = 0; c < d1; ++c) col[r + c * d0] = row[r * d1 + c];
  } else {
    long long j; for (j = 0; j < n; ++j) col[j] = row[j];
  }
}

/* Byte size of one element of ONNX type et; 0 if et is non-numeric (unsupported). */
static size_t casadi_onnxruntime_elem_size(long long et) {
  switch (et) {
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_DOUBLE:
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT64:
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT64:  return 8;
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT:
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT32:
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT32:  return 4;
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT16:
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT16:
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT16:
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_BFLOAT16: return 2;
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT8:
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT8:
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_BOOL:    return 1;
    default: return 0;
  }
}

/* Pack a row-major double buffer into the caller-owned tensor buffer dst of ONNX type et.
 * dst must hold casadi_onnxruntime_elem_size(et)*n bytes. Returns 0 on success, 1 if non-numeric. */
#define CASADI_ORT_PACK(CTYPE) \
  { CTYPE* b = (CTYPE*) dst; for (k = 0; k < n; ++k) b[k] = (CTYPE) row[k]; return 0; }
static int casadi_onnxruntime_pack_into(long long et, const double* row, long long n, void* dst) {
  long long k;
  switch (et) {
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT:  CASADI_ORT_PACK(float)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_DOUBLE: CASADI_ORT_PACK(double)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT8:   CASADI_ORT_PACK(int8_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT8:  CASADI_ORT_PACK(uint8_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT16:  CASADI_ORT_PACK(int16_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT16: CASADI_ORT_PACK(uint16_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT32:  CASADI_ORT_PACK(int32_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT32: CASADI_ORT_PACK(uint32_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT64:  CASADI_ORT_PACK(int64_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT64: CASADI_ORT_PACK(uint64_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_BOOL:   CASADI_ORT_PACK(uint8_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT16: {
      uint16_t* b = (uint16_t*) dst;
      for (k = 0; k < n; ++k) b[k] = casadi_onnxruntime_d2half(row[k]);
      return 0; }
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_BFLOAT16: {
      uint16_t* b = (uint16_t*) dst;
      for (k = 0; k < n; ++k) b[k] = casadi_onnxruntime_d2bf16(row[k]);
      return 0; }
    default: return 1;
  }
}
#undef CASADI_ORT_PACK

/* Unpack an ONNX tensor of type et into a row-major double buffer. Returns 0 on success. */
#define CASADI_ORT_UNPACK(CTYPE) \
  { const CTYPE* b = (const CTYPE*) td; for (k = 0; k < n; ++k) row[k] = (double) b[k]; return 0; }
static int casadi_onnxruntime_unpack(long long et, const void* td, long long n, double* row) {
  long long k;
  switch (et) {
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT:  CASADI_ORT_UNPACK(float)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_DOUBLE: CASADI_ORT_UNPACK(double)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT8:   CASADI_ORT_UNPACK(int8_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT8:  CASADI_ORT_UNPACK(uint8_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT16:  CASADI_ORT_UNPACK(int16_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT16: CASADI_ORT_UNPACK(uint16_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT32:  CASADI_ORT_UNPACK(int32_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT32: CASADI_ORT_UNPACK(uint32_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT64:  CASADI_ORT_UNPACK(int64_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT64: CASADI_ORT_UNPACK(uint64_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_BOOL:   CASADI_ORT_UNPACK(uint8_t)
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT16:
      { const uint16_t* b = (const uint16_t*) td;
        for (k = 0; k < n; ++k) row[k] = casadi_onnxruntime_half2d(b[k]);
        return 0; }
    case ONNX_TENSOR_ELEMENT_DATA_TYPE_BFLOAT16:
      { const uint16_t* b = (const uint16_t*) td;
        for (k = 0; k < n; ++k) row[k] = casadi_onnxruntime_bf162d(b[k]);
        return 0; }
    default: return 1;
  }
}
#undef CASADI_ORT_UNPACK

/* Create the session (idempotent). Returns 0 on success. */
static int casadi_onnxruntime_init(struct casadi_onnxruntime_data* d,
                                   const struct casadi_onnxruntime_prob* p) {
  const OrtApi* api = casadi_onnxruntime_api();
  OrtSessionOptions* so = 0;
  if (d->session) return 0;
  if (!api) return 1;
  if (api->CreateEnv(ORT_LOGGING_LEVEL_WARNING, "casadi", &d->env)) return 1;
  if (api->CreateSessionOptions(&so)) return 1;
  if (api->CreateSessionFromArray(d->env, p->model_data, (size_t)p->model_size, so, &d->session)) {
    api->ReleaseSessionOptions(so);
    return 1;
  }
  api->ReleaseSessionOptions(so);
  if (api->CreateCpuMemoryInfo(OrtArenaAllocator, OrtMemTypeDefault, &d->mem)) return 1;
  return 0;
}

/* Build the value-independent per-eval scaffolding once (idempotent): the typed input buffers,
 * a reusable input tensor (OrtValue) viewing each buffer, the conversion scratch, and the
 * constant (baked/unwired) input values. Requires the session/mem from casadi_onnxruntime_init.
 * This is the per-instance half of the ONNX "prepare" step; solve() then does no allocation.
 * Returns 0 on success. */
static int casadi_onnxruntime_prepare(struct casadi_onnxruntime_data* d,
                                      const struct casadi_onnxruntime_prob* p) {
  const OrtApi* api = casadi_onnxruntime_api();
  long long i, off, voff, k, maxnel = 1;
  if (d->prepared) return 0;
  if (!api) return 1;
  /* one scratch row, sized to the largest input or output tensor */
  for (i = 0; i < p->n_in;  ++i) if (p->in_numel[i]  > maxnel) maxnel = p->in_numel[i];
  for (i = 0; i < p->n_out; ++i) if (p->out_numel[i] > maxnel) maxnel = p->out_numel[i];
  d->row  = (double*)    malloc(sizeof(double)   * (size_t)maxnel);
  d->inv  = (OrtValue**) calloc((size_t)p->n_in,  sizeof(OrtValue*));
  d->outv = (OrtValue**) calloc((size_t)p->n_out, sizeof(OrtValue*));
  d->buf  = (void**)     calloc((size_t)p->n_in,  sizeof(void*));
  if (!d->row || !d->inv || !d->outv || !d->buf) return 1;
  for (i = 0, off = 0, voff = 0; i < p->n_in;
       off += p->in_ndim[i], voff += p->in_numel[i], ++i) {
    long long nel = p->in_numel[i], nd = p->in_ndim[i], src = p->in_src[i];
    size_t esz = casadi_onnxruntime_elem_size(p->in_elem_type[i]);
    int64_t* shp;
    if (esz == 0) return 1;
    d->buf[i] = malloc(esz * (size_t)nel);
    if (!d->buf[i]) return 1;
    /* Inputs whose values never change are packed once now; arg-fed inputs are filled per solve() */
    if (src == -2) {  /* baked value */
      casadi_onnxruntime_to_row(nd, p->in_dims + off, p->in_val + voff, d->row, nel);
      casadi_onnxruntime_pack_into(p->in_elem_type[i], d->row, nel, d->buf[i]);
    } else if (src < 0) {  /* unwired -> NaN (poisons dependent outputs) where the type allows, else 0 */
      double fill = casadi_onnxruntime_is_float(p->in_elem_type[i]) ? (double) NAN : 0.0;
      for (k = 0; k < nel; ++k) d->row[k] = fill;
      casadi_onnxruntime_pack_into(p->in_elem_type[i], d->row, nel, d->buf[i]);
    }
    shp = (int64_t*) malloc(sizeof(int64_t) * (size_t)(nd > 0 ? nd : 1));
    for (k = 0; k < nd; ++k) shp[k] = (int64_t) p->in_dims[off + k];
    if (api->CreateTensorWithDataAsOrtValue(d->mem, d->buf[i], esz * (size_t)nel, shp, (size_t)nd,
          (ONNXTensorElementDataType) p->in_elem_type[i], &d->inv[i])) { free(shp); return 1; }
    free(shp);
  }
  d->prepared = 1;
  return 0;
}

/* Evaluate the model. arg holds the exposed inputs; res the exposed outputs. Returns 0 on success.
 * Allocation and tensor creation happen once in casadi_onnxruntime_prepare; per call this only
 * refills the arg-fed input buffers, runs the pre-created input tensors, and unpacks the outputs. */
static int casadi_onnxruntime_solve(struct casadi_onnxruntime_data* d,
                                    const struct casadi_onnxruntime_prob* p,
                                    const double** arg, double** res) {
  const OrtApi* api = casadi_onnxruntime_api();
  long long i, off, k, ret = 1;
  if (!api || casadi_onnxruntime_prepare(d, p)) return 1;

  /* Refill only the inputs fed by the caller; baked/unwired buffers keep their prepared values */
  for (i = 0, off = 0; i < p->n_in; off += p->in_ndim[i], ++i) {
    long long nel = p->in_numel[i], nd = p->in_ndim[i], src = p->in_src[i];
    if (src < 0) continue;
    if (arg[src]) {
      casadi_onnxruntime_to_row(nd, p->in_dims + off, arg[src], d->row, nel);
    } else {
      /* declared input not supplied -> NaN/0 sentinel, matching the unwired convention */
      double fill = casadi_onnxruntime_is_float(p->in_elem_type[i]) ? (double) NAN : 0.0;
      for (k = 0; k < nel; ++k) d->row[k] = fill;
    }
    casadi_onnxruntime_pack_into(p->in_elem_type[i], d->row, nel, d->buf[i]);
  }

  if (api->Run(d->session, 0, p->input_names, (const OrtValue* const*) d->inv, (size_t)p->n_in,
               p->output_names, (size_t)p->n_out, d->outv)) goto cleanup;

  for (i = 0, off = 0; i < p->n_out; off += p->out_ndim[i], ++i) {
    void* td;
    if (res[i]) {
      if (api->GetTensorMutableData(d->outv[i], &td)) goto cleanup;
      if (casadi_onnxruntime_unpack(p->out_elem_type[i], td, p->out_numel[i], d->row)) goto cleanup;
      casadi_onnxruntime_from_row(p->out_ndim[i], p->out_dims + off, d->row, res[i], p->out_numel[i]);
    }
    api->ReleaseValue(d->outv[i]); d->outv[i] = 0;  /* Run owns these; release after reading */
  }
  ret = 0;

cleanup:
  for (i = 0; i < p->n_out; ++i) if (d->outv[i]) { api->ReleaseValue(d->outv[i]); d->outv[i] = 0; }
  return (int) ret;
}

#endif /* CASADI_ORT_RUNTIME_H */
