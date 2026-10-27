package org.casadi;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.File;
import java.net.URL;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Stream;

/**
 * Loads the casadi_java native library, from the first of:
 * <ol>
 * <li>the directory in system property {@code casadi.library.path} or environment variable {@code CASADI_JAVA_LIBRARY_PATH};</li>
 * <li>a natives jar on the class path ({@code org/casadi/natives/<platform>/}), extracted once to a cache directory
 *     (system property {@code casadi.cache.dir}, default under the user cache directory);</li>
 * <li>the directory of casadi.jar (an installed casadi);</li>
 * <li>{@code java.library.path}.</li>
 * </ol>
 * Plugins are looked up next to the library that was loaded.
 */
final class NativeLoader {
  /** Directory of the loaded casadi_java library, null if unknown */
  static String directory;

  /** casadi headers for code generation: include next to the library, or ../include */
  static String includeDirectory() {
    if (directory == null) return null;
    Path d = Paths.get(directory);
    Path inc = d.resolve("include");
    if (!Files.isDirectory(inc) && d.getParent() != null) inc = d.getParent().resolve("include");
    return inc.toString();
  }

  private NativeLoader() {}

  static synchronized void load() {
    if (directory != null) return;
    String lib = System.mapLibraryName("casadi_java");
    String dir = System.getProperty("casadi.library.path", System.getenv("CASADI_JAVA_LIBRARY_PATH"));
    if (dir != null && !dir.isEmpty()) {
      loadFrom(Paths.get(dir), lib);
      return;
    }
    String res = "natives/" + platform() + "/";
    URL index = NativeLoader.class.getResource(res + "files.txt");
    if (index != null) {
      loadFrom(extract(res, index), lib);
      return;
    }
    Path here = codeDirectory();
    if (here != null && Files.isRegularFile(here.resolve(lib))) {
      loadFrom(here, lib);
      return;
    }
    System.loadLibrary("casadi_java");
    for (String p : System.getProperty("java.library.path", "").split(File.pathSeparator)) {
      if (!p.isEmpty() && Files.isRegularFile(Paths.get(p, lib))) {
        directory = Paths.get(p).toAbsolutePath().toString();
        break;
      }
    }
  }

  /** e.g. linux-x86_64, linux-aarch64, windows-x86_64, macos-arm64, macos-x86_64 */
  static String platform() {
    String os = System.getProperty("os.name").toLowerCase();
    String arch = System.getProperty("os.arch").toLowerCase();
    os = os.startsWith("windows") ? "windows" : os.startsWith("mac") ? "macos" : os.replace(' ', '_');
    if (arch.equals("amd64") || arch.equals("x86_64")) {
      arch = "x86_64";
    } else if (arch.equals("aarch64") || arch.equals("arm64")) {
      arch = os.equals("macos") ? "arm64" : "aarch64";
    }
    return os + "-" + arch;
  }

  private static boolean isWindows() {
    return System.getProperty("os.name").toLowerCase().startsWith("windows");
  }

  private static void loadFrom(Path dir, String lib) {
    dir = dir.toAbsolutePath();
    Path p = dir.resolve(lib);
    if (!Files.isRegularFile(p)) {
      throw new UnsatisfiedLinkError("casadi: " + p + " not found");
    }
    // Windows binds dependencies by module name, searching PATH rather than dir: preload ours, leaves first
    if (isWindows()) preloadImports(dir, lib, new HashSet<>());
    System.load(p.toString());
    directory = dir.toString();
  }

  /** Loads the DLLs in dir that lib imports, transitively, dependencies before dependents */
  private static void preloadImports(Path dir, String lib, Set<String> seen) {
    for (String dep : imports(dir.resolve(lib))) {
      Path d = dir.resolve(dep);
      if (!seen.add(dep.toLowerCase()) || !Files.isRegularFile(d)) continue;
      preloadImports(dir, dep, seen);
      System.load(d.toString());
    }
  }

  /** DLL names in the import table of a PE file, empty if unreadable */
  static List<String> imports(Path pe) {
    List<String> r = new ArrayList<>();
    try {
      ByteBuffer b = ByteBuffer.wrap(Files.readAllBytes(pe)).order(ByteOrder.LITTLE_ENDIAN);
      int nt = b.getInt(0x3c);
      if (b.getInt(nt) != 0x4550) return r;
      int nsec = b.getShort(nt + 6) & 0xffff;
      int opt = nt + 24, optSize = b.getShort(nt + 20) & 0xffff;
      int dirs = opt + ((b.getShort(opt) & 0xffff) == 0x20b ? 112 : 96);
      int importRva = b.getInt(dirs + 8);
      if (importRva == 0) return r;
      int sec = opt + optSize;
      for (int d = rvaToOffset(b, sec, nsec, importRva); b.getInt(d + 12) != 0; d += 20) {
        int o = rvaToOffset(b, sec, nsec, b.getInt(d + 12));
        StringBuilder n = new StringBuilder();
        while (b.get(o) != 0) n.append((char) b.get(o++));
        r.add(n.toString());
      }
    } catch (IOException | RuntimeException e) {
      // not a readable PE file: nothing to preload
    }
    return r;
  }

  private static int rvaToOffset(ByteBuffer b, int sec, int nsec, int rva) {
    for (int i = 0; i < nsec; i++, sec += 40) {
      int va = b.getInt(sec + 12), size = Math.max(b.getInt(sec + 8), b.getInt(sec + 16));
      if (rva >= va && rva < va + size) return rva - va + b.getInt(sec + 20);
    }
    throw new IllegalArgumentException("rva outside sections");
  }

  private static Path codeDirectory() {
    try {
      Path p = Paths.get(NativeLoader.class.getProtectionDomain().getCodeSource().getLocation().toURI());
      return Files.isDirectory(p) ? p : p.getParent();
    } catch (Exception e) {
      return null;
    }
  }

  private static Path cacheRoot() {
    String d = System.getProperty("casadi.cache.dir");
    if (d != null && !d.isEmpty()) return Paths.get(d);
    String home = System.getProperty("user.home");
    if (isWindows()) {
      String a = System.getenv("LOCALAPPDATA");
      return Paths.get(a != null ? a : home, "casadi-java");
    }
    if (System.getProperty("os.name").toLowerCase().startsWith("mac")) {
      return Paths.get(home, "Library", "Caches", "casadi-java");
    }
    String x = System.getenv("XDG_CACHE_HOME");
    return x != null && !x.isEmpty() ? Paths.get(x, "casadi-java") : Paths.get(home, ".cache", "casadi-java");
  }

  /** Extracts the natives listed in index, once per index content */
  private static Path extract(String res, URL index) {
    try {
      byte[] idx;
      try (InputStream in = index.openStream()) {
        idx = in.readAllBytes();
      }
      Path root = cacheRoot();
      Path target = root.resolve(hex(idx).substring(0, 16));
      if (Files.isRegularFile(target.resolve(".complete"))) return target;
      Files.createDirectories(root);
      Path tmp = Files.createTempDirectory(root, ".extract");
      try {
        BufferedReader r = new BufferedReader(new InputStreamReader(new java.io.ByteArrayInputStream(idx), StandardCharsets.UTF_8));
        for (String line; (line = r.readLine()) != null;) {
          String[] f = line.trim().split(" ");
          if (f.length == 3 && f[0].equals("file")) {
            try (InputStream in = NativeLoader.class.getResourceAsStream(res + f[2])) {
              if (in == null) throw new IOException("missing resource " + res + f[2]);
              Files.copy(in, tmp.resolve(f[2]));
            }
          }
        }
        r = new BufferedReader(new InputStreamReader(new java.io.ByteArrayInputStream(idx), StandardCharsets.UTF_8));
        for (String line; (line = r.readLine()) != null;) {
          String[] f = line.trim().split(" ");
          if (f.length == 3 && f[0].equals("link")) {
            try {
              Files.createSymbolicLink(tmp.resolve(f[1]), Paths.get(f[2]));
            } catch (IOException | UnsupportedOperationException e) {
              Files.copy(tmp.resolve(f[2]), tmp.resolve(f[1]));
            }
          }
        }
        Files.createFile(tmp.resolve(".complete"));
        try {
          Files.move(tmp, target, StandardCopyOption.ATOMIC_MOVE);
        } catch (FileAlreadyExistsException e) {
          // extracted concurrently
        } catch (IOException e) {
          if (!Files.isRegularFile(target.resolve(".complete"))) throw e;
        }
      } finally {
        if (Files.exists(tmp)) delete(tmp);
      }
      return target;
    } catch (IOException e) {
      throw new UnsatisfiedLinkError("casadi: cannot extract native libraries: " + e);
    }
  }

  private static void delete(Path p) {
    try (Stream<Path> s = Files.walk(p)) {
      s.sorted(Comparator.reverseOrder()).forEach(f -> f.toFile().delete());
    } catch (IOException e) {
      // best effort
    }
  }

  private static String hex(byte[] b) {
    try {
      StringBuilder s = new StringBuilder();
      for (byte x : MessageDigest.getInstance("SHA-256").digest(b)) s.append(String.format("%02x", x));
      return s.toString();
    } catch (NoSuchAlgorithmException e) {
      throw new IllegalStateException(e);
    }
  }
}
