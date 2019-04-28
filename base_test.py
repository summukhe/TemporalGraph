from temporal_graph import *

if __name__ == "__main__":
    edges = [('A', 'B'), ('A', 'E'), ('A', 'C'),
             ('D', 'B'), ('E', 'D'), ('B', 'E'),
             ('G', 'H'), ('H', 'I'), ('G', 'K')]
    g = WeightedGraph(directed=False)

    print(dir(g))
    for v1, v2 in edges:
        g.add_edge(v1, v2)
    print(g.adjacency())

    pulse = PulseSignal(ts=0, te=4)
    process = InformationDiffusionOnStaticGraph()
    process.start(g, start_nodes=['A', 'C'], pulse=pulse)
    for i in range(30):
        process.step()
    print(process.get_history('H'))