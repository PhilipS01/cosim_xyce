// ELK (Eclipse Layout Kernel via elkjs) layout for the schematic renderer.
// stdin  : {"components":[{"id","a","b"}...], "nets":[net...], "port":..., "gnd":...}
//          each component is a 2-terminal element between electrical nets a and b.
// stdout : {"nets":{net:[x,y]...}, "routes":{compId:[[x,y]...]...}, "size":[w,h]}
//          nets = placed net (junction) positions; routes = the orthogonal wire polyline ELK routed
//          for each component edge (start .. bendpoints .. end).
const ELK = require('elkjs');
const elk = new ELK();
let input = '';
process.stdin.on('data', d => input += d);
process.stdin.on('end', async () => {
  const g = JSON.parse(input);
  const nodeOpts = id => {
    const o = {};
    if (id === g.port) o['elk.layered.layering.layerConstraint'] = 'FIRST_SEPARATE';
    if (id === g.gnd)  o['elk.layered.layering.layerConstraint'] = 'LAST_SEPARATE';
    return o;
  };
  const graph = {
    id: 'root',
    layoutOptions: {
      'elk.algorithm': 'layered',
      'elk.direction': 'DOWN',
      'elk.edgeRouting': 'ORTHOGONAL',
      'elk.layered.spacing.nodeNodeBetweenLayers': '70',
      'elk.spacing.nodeNode': '75',
      'elk.spacing.edgeNode': '30',
      'elk.spacing.edgeEdge': '20',
      'elk.layered.spacing.edgeNodeBetweenLayers': '30',
      'elk.layered.nodePlacement.strategy': 'NETWORK_SIMPLEX',
    },
    // Net nodes are ~point-sized so every wire on a net converges to the SAME coordinate -> clean
    // orthogonal junctions (a 12x12 box would let routes hit different border points -> diagonal stubs).
    children: g.nets.map(id => ({ id, width: 1, height: 1, layoutOptions: nodeOpts(id) })),
    edges: g.components.map(c => ({ id: c.id, sources: [c.a], targets: [c.b] })),
  };
  const res = await elk.layout(graph);
  const nets = {};
  res.children.forEach(n => nets[n.id] = [n.x + n.width / 2, n.y + n.height / 2]);
  const routes = {};
  (res.edges || []).forEach(e => {
    const s = (e.sections || [])[0];
    if (!s) { routes[e.id] = []; return; }
    const pts = [[s.startPoint.x, s.startPoint.y]];
    (s.bendPoints || []).forEach(b => pts.push([b.x, b.y]));
    pts.push([s.endPoint.x, s.endPoint.y]);
    routes[e.id] = pts;
  });
  process.stdout.write(JSON.stringify({ nets, routes, size: [res.width, res.height] }));
});
