package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.casadi.DM;
import org.casadi.DeserializerBase;
import org.casadi.Function;
import org.casadi.SerializerBase;
import org.casadi.StringDeserializer;
import org.casadi.StringSerializer;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** Port of test/python/serialize.py */
public class SerializeTest extends CasadiTestCase {


  @Test
  public void test_compat() {
    File dir = new File("../python/serialize_3.5.5");
    Assumptions.assumeTrue(dir.isDirectory(), "no serialize_3.5.5 corpus");
    Map<String, List<String>> inputs = new TreeMap<>(), outputs = new TreeMap<>();
    String[] names = dir.list();
    List<String> functions = new ArrayList<>();
    for (String n : names) {
      if (n.endsWith(".casadi")) functions.add(n);
      (n.contains("_in") ? inputs : outputs).computeIfAbsent(n.split("_")[0], k -> new ArrayList<>()).add(n);
    }
    for (String fun : functions) {
      Function f;
      try {
        f = Function.load(new File(dir, fun).getPath());
      } catch (RuntimeException e) {
        String msg = String.valueOf(e.getMessage());
        if (msg.contains("DllLibrary::init_handle")) continue;
        if (msg.contains("DeSerialization of Integrator failed")) continue;
        if (msg.contains("DeSerialization of Rootfinder failed")) continue;
        boolean skip = false;
        for (String token : msg.split("\\s+")) {
          if (token.startsWith("libcasadi_")) {
            String[] parts = token.replaceAll("^['\"]+|['\"]+$", "").split("_", 3);
            if (parts.length == 3 && !hasPlugin(parts[1], parts[2])) {
              skip = true;
              break;
            }
          }
        }
        if (skip) continue;
        throw e;
      }
      String key = fun.split("\\.")[0];
      List<String> in = new ArrayList<>(inputs.getOrDefault(key, List.of()));
      List<String> out = new ArrayList<>(outputs.getOrDefault(key, List.of()));
      in.sort(null);
      out.sort(null);
      for (int k = 0; k < Math.min(in.size(), out.size()); k++) {
        List<DM> inp = f.generate_in(new File(dir, in.get(k)).getPath());
        List<DM> outp_ref = f.generate_out(new File(dir, out.get(k)).getPath());
        List<DM> outp = f.call(inp);
        assertEquals(outp_ref.size(), outp.size());
        for (int i = 0; i < outp.size(); i++) {
          double tol = 1e-15;
          String name = f.name_out(i);
          if (f.toString().contains("SuperscsInterface")) {
            if (name.startsWith("lam_")) continue;
            tol = 1e-7;
          }
          String os = System.getProperty("os.name").toLowerCase();
          if (os.contains("mac") || os.contains("win")) tol = 1e-7;
          checkarray(outp.get(i), outp_ref.get(i), tol, fun + "/" + f);
        }
      }
    }
  }

  private static boolean hasPlugin(String kind, String plugin) {
    switch (kind) {
      case "nlpsol": return hasNlpsol(plugin);
      case "conic": return hasConic(plugin);
      case "integrator": return hasIntegrator(plugin);
      case "rootfinder": return hasRootfinder(plugin);
      case "linsol": return hasLinsol(plugin);
      case "interpolant": return org.casadi.ca.has_interpolant(plugin);
      default: return true;
    }
  }

  private static void checkEqual(Object a, Object b) {
    if (a instanceof Map) {
      Map<?, ?> x = (Map<?, ?>) a, y = (Map<?, ?>) b;
      assertEquals(new TreeMap<>(x).keySet(), new TreeMap<>(y).keySet());
      for (Object k : x.keySet()) checkEqual(x.get(k), y.get(k));
    } else if (a instanceof List) {
      List<?> x = (List<?>) a, y = (List<?>) b;
      assertEquals(x.size(), y.size());
      for (int i = 0; i < x.size(); i++) checkEqual(x.get(i), y.get(i));
    } else {
      assertEquals(a, b);
    }
  }

  @Test
  public void test_identity() {
    List<Object> obj = Arrays.asList("foo", Map.of("foo", "bar"), List.of(Map.of("a", 3L), Map.of("b", 9L)),
        List.of("a", 5L), Map.of("foo", List.of("a", 5L)), Map.of("foo", List.of(List.of("a", 5L), List.of("b", 2L))),
        List.of(List.of("a", 5L), List.of("b", 2L)), List.of(List.of(List.of("a", 5L), List.of("b", 2L))));
    for (Object e : obj) {
      StringSerializer ss = new StringSerializer();
      if (e instanceof String) ss.pack((String) e);
      else ss.pack(e);
      StringDeserializer ds = new StringDeserializer(ss.encode());
      Object r = ds.unpack();
      System.out.println(e + " " + r);
      checkEqual(e, r);
    }
  }

  // test_pickling_context: python-only (pickle protocol, global_pickle_context)
}
