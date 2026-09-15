// Asynchronous plugin entry points must be usable as promises.
import createCasadi from 'casadi';
async function checkLoaderTypes() {
  const ca = await createCasadi();
  const pending: Promise<void>[] = [
    ca.load_nlpsol('fatrop'), ca.load_conic('osqp'), ca.load_blas('classic'),
    ca.load_filesystem('ghc'), ca.Linsol.load_plugin('csparsecholesky'),
    ca.XmlFile.load_plugin('tinyxml'), ca.Importer.load_plugin('shell'),
    ca.ModelicaParser.load_plugin('lacemodelica'),
  ];
  await Promise.all(pending);
}
