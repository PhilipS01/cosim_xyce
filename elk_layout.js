const ELK = require('elkjs');
const elk = new ELK();
let input = '';
process.stdin.on('data', d => input += d);
process.stdin.on('end', async () => {
  const g = JSON.parse(input);                    // {nodes, edges, port, gnd}
  const opt = id => {
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
      'elk.layered.spacing.nodeNodeBetweenLayers': '40',
      'elk.spacing.nodeNode': '40',
    },
    children: g.nodes.map(id => ({ id, width: 20, height: 20, layoutOptions: opt(id) })),
    edges: g.edges.map((e, i) => ({ id: 'e' + i, sources: [e[0]], targets: [e[1]] })),
  };
  const res = await elk.layout(graph);
  const pos = {};
  res.children.forEach(c => pos[c.id] = [c.x, c.y]);
  process.stdout.write(JSON.stringify(pos));
});
