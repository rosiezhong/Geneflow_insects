from pathlib import Path
import tskit
import sys
import numpy as np

dir = Path(sys.argv[1])
num_samples = int(sys.argv[2])
out = sys.argv[3]


def sum_tree_span(ts, tree_dict):
    tracked0 = tskit.Tree(ts, tracked_samples=[0, 1])
    tracked1 = tskit.Tree(ts, tracked_samples=[2, 3])
    tracked2 = tskit.Tree(ts, tracked_samples=[0, 2])
    tracked3 = tskit.Tree(ts, tracked_samples=[1, 3])
    tracked4 = tskit.Tree(ts, tracked_samples=[0, 3])
    tracked5 = tskit.Tree(ts, tracked_samples=[1, 2])

    for tree_index in range(ts.num_trees):

        tracked0.next()
        tracked1.next()
        tracked2.next()
        tracked3.next()
        tracked4.next()
        tracked5.next()
        assert tree_index == tracked0.index
        assert tree_index == tracked1.index
        assert tree_index == tracked2.index
        assert tree_index == tracked3.index
        assert tree_index == tracked4.index
        assert tree_index == tracked5.index

        if any(
            tracked0.num_children(u) != 2
            for u in tracked0.nodes()
            if tracked0.is_internal(u)
        ):  # Ignore trees with polytomies
            continue
        if not tracked0.has_single_root:  # Ignore trees with multiple roots
            continue

        NODES = list(tracked0.nodes(order="timeasc"))
        assert NODES[-1] == tracked0.root
        node4, node5 = NODES[4:6]

        counts0 = [
            tracked0.num_tracked_samples(node4),
            tracked0.num_tracked_samples(node5),
        ]
        counts1 = [
            tracked1.num_tracked_samples(node4),
            tracked1.num_tracked_samples(node5),
        ]
        counts2 = [
            tracked2.num_tracked_samples(node4),
            tracked2.num_tracked_samples(node5),
        ]
        counts3 = [
            tracked3.num_tracked_samples(node4),
            tracked3.num_tracked_samples(node5),
        ]
        counts4 = [
            tracked4.num_tracked_samples(node4),
            tracked4.num_tracked_samples(node5),
        ]
        counts5 = [
            tracked5.num_tracked_samples(node4),
            tracked5.num_tracked_samples(node5),
        ]

        counts0_set = set(counts0)
        counts1_set = set(counts1)
        counts2_set = set(counts2)
        counts3_set = set(counts3)
        counts4_set = set(counts4)
        counts5_set = set(counts5)

        # Span of 18 distinct dtwfescent trees
        # BALANCED
        if counts0_set == counts1_set == {0, 2}:  # b0123
            if tracked0.time(node4) == tracked0.time(node5):
                tree_dict["b0123_eq"] += tracked0.span
            else:
                if counts0[0] == 2:
                    tree_dict["b0123_01"] += tracked0.span
                if counts0[0] == 0:
                    tree_dict["b0123_23"] += tracked0.span
        if counts2_set == counts3_set == {0, 2}:  # b0213
            if tracked0.time(node4) == tracked0.time(node5):
                tree_dict["b0213_eq"] += tracked0.span
            else:
                if counts2[0] == 2:
                    tree_dict["b0213_02"] += tracked0.span
                if counts3[0] == 0:
                    tree_dict["b0213_13"] += tracked0.span
        if counts4_set == counts5_set == {0, 2}:  # b0312
            if tracked0.time(node4) == tracked0.time(node5):
                tree_dict["b0312_eq"] += tracked0.span
            else:
                if counts4[0] == 2:
                    tree_dict["b0312_03"] += tracked0.span
                if counts5[0] == 0:
                    tree_dict["b0312_12"] += tracked0.span
        # UNBALANCED
        if counts0_set == {2}:  # u01XX
            if counts2_set == counts5_set == {1, 2}:
                tree_dict["u0123"] += tracked0.span
            if counts3_set == counts4_set == {1, 2}:
                tree_dict["u0132"] += tracked0.span
        if counts1_set == {2}:  # u23XX
            if counts2_set == counts4_set == {1, 2}:
                tree_dict["u2301"] += tracked0.span
            if counts3_set == counts5_set == {1, 2}:
                tree_dict["u2310"] += tracked0.span
        if counts2_set == {2}:  # u02XX
            if counts0_set == counts5_set == {1, 2}:
                tree_dict["u0213"] += tracked0.span
            if counts1_set == counts4_set == {1, 2}:
                tree_dict["u0231"] += tracked0.span
        if counts3_set == {2}:  # u13XX
            if counts0_set == counts4_set == {1, 2}:
                tree_dict["u1302"] += tracked0.span
            if counts1_set == counts5_set == {1, 2}:
                tree_dict["u1320"] += tracked0.span
        if counts4_set == {2}:  # u03XX
            if counts0_set == counts3_set == {1, 2}:
                tree_dict["u0312"] += tracked0.span
            if counts1_set == counts2_set == {1, 2}:
                tree_dict["u0321"] += tracked0.span
        if counts5_set == {2}:  # u12XX
            if counts0_set == counts2_set == {1, 2}:
                tree_dict["u1203"] += tracked0.span
            if counts1_set == counts3_set == {1, 2}:
                tree_dict["u1230"] += tracked0.span

    return tree_dict


def get_tree_dict():
    tree_dict = {
        "b0213_eq": 0,
        "b0312_eq": 0,
        "b0123_eq": 0,
        "b0213_02": 0,
        "b0213_13": 0,
        "b0312_03": 0,
        "b0312_12": 0,
        "u0123": 0,
        "u0132": 0,
        "u2301": 0,
        "u2310": 0,
        "b0123_01": 0,
        "b0123_23": 0,
        "u0213": 0,
        "u0312": 0,
        "u1203": 0,
        "u1302": 0,
        "u0231": 0,
        "u0321": 0,
        "u1230": 0,
        "u1320": 0,
    }
    return tree_dict


def get_roman_numerals():
    roman_numerals = [
        "I",
        "II",
        "III",
        "IV",
        "V",
        "VI",
        "VII",
        "VIII",
        "IX",
        "X",
        "XI",
        "XII",
        "XIII",
        "XIV",
        "XV",
        "XVI",
        "XVII",
        "XVIII",
        "XIX",
        "XX",
        "XXI",
    ]
    return roman_numerals


def subsample_ts_by_nodes(ts, ind_num):
    nodes1 = np.random.choice(
        [n for n in range(ts.num_samples) if n % 2 == 0], size=ind_num, replace=False
    )
    nodes2 = nodes1 + 1
    nodes = np.concatenate([nodes1, nodes2])
    nodes = np.sort(nodes)
    ts = ts.simplify(nodes)
    return ts


tree_dict = get_tree_dict()
for ts in dir.glob("*.trees"):
    ts = tskit.load(str(ts))
    ts = subsample_ts_by_nodes(ts, num_samples)
    tree_dict = sum_tree_span(ts, tree_dict)

roman_numerals = get_roman_numerals()
tree_dict = dict(zip(roman_numerals, tree_dict.values()))


with open(out, "w") as f:
    for key, value in tree_dict.items():
        f.write(f"{key}\t{value}\n")
