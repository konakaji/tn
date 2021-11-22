def _four_loop(self, network, all_pairs):
    indices = set()  # the indices of pairs that is not one-loop
    n_loop = 0
    for i, pair in enumerate(all_pairs):
        if PlugUtil.connected(pair[0], pair[1]):
            network.remove_edge(pair[0])
            network.remove_edge(pair[1])
            n_loop = n_loop + 1
            indices.add(i)
    return indices, n_loop


map = {1: Factor.D, 2: Factor.D2, 3: Factor.D3, 4: Factor.D4}
factors = [Factor.G]
factor_count = len(all_pairs) - len(pairs)
if factor_count != 0:
    factors.append(map[factor_count])
for pair in pairs:
    l_plug, r_plug = pair
    lp = l_plug.edge.left_plug
    rp = r_plug.edge.right_plug
    network.remove_edge(l_plug)
    network.remove_edge(r_plug)
    network.add_edge(lp, rp)


def _two_loop(self, network, all_pairs):
    for i in [0, 1]:
        if not (PlugUtil.connected(all_pairs[i][0], all_pairs[2 + i][1])
                and PlugUtil.connected(all_pairs[i][1], all_pairs[2 + i][0])):
            return {}, 0
    indices = set()
    for i in [0, 1]:
        network.remove_edge(all_pairs[i][0])
        network.remove_edge(all_pairs[2 + i][1])
        network.remove_edge(all_pairs[i][1])
        network.remove_edge(all_pairs[2 + i][0])
        indices.add(i)
        indices.add(2 + i)
    return indices, 1

# def _four_loop(self, network, all_pairs):
#     indices = set()  # the indices of pairs that is not one-loop
#     n_loop = 0
#     if PlugUtil.connected(all_pairs[0][0], all_pairs[2 + i][1]) \
#             and PlugUtil.connected(all_pairs[i][1], all_pairs[2 + i][0]):
#         network.remove_edge(pair[0])
#         network.remove_edge(pair[1])
#         n_loop = n_loop + 1
#         indices.add(i)
#     return indices, n_loop

    for n1, n2 in zip(sorted(ns[0].node_map.items()), sorted(ns[1].node_map.items())):
        if n1.__hash__() != n2.__hash__():
            print("hoge")
    for i1, i2 in zip(sorted(ns[0].edge_map.items()), sorted(ns[1].edge_map.items())):
        n = i1[1]
        n2 = i2[1]
        if n.__hash__() != n2.__hash__():
            print(n.__hash__(), n2.__hash__())
            e: Edge = n
            e2: Edge = n2
            if e.left_plug.__hash__() != e2.left_plug.__hash__():
                print("left")
                print(ns[0].node_map[ns[0].node_map[e.left_plug.node_id].id]._location)
                print(ns[1].node_map[ns[1].node_map[e2.left_plug.node_id].id]._location)
            if e.right_plug.node_id.__hash__() != e2.right_plug.node_id.__hash__():
                print("right")
                print(ns[0].node_map[ns[0].node_map[e.right_plug.node_id].id]._location)
                print(ns[1].node_map[ns[1].node_map[e2.right_plug.node_id].id]._location)
