
from dataclasses import dataclass
from pprint import pprint
from collections import Counter
from typing import Optional

import numpy

import toytree
import toyplot
import toyplot.svg
import toyplot.png


@dataclass
class MigrationEvent:
    origin_node_id: int
    destination_node_id: int
    weight: Optional[float] = None
    support: Optional[int] = None


class TreeWithMigrations:
    def __init__(self, tree, migration_events=None):
        self.tree = tree
        if migration_events is None:
            migration_events = []
        self.migration_events = migration_events


class TreeWithMigrationsList:
    def __init__(self, trees_with_migrations=None):
        if trees_with_migrations is None:
            trees_with_migrations = []

        self.trees = trees_with_migrations

    def append(self, tree):
        self.trees.append(tree)

    def extend(self, trees):
        self.trees.extend(trees)

    def get_consensus_tree(self):
        only_trees = toytree.mtree([tree.tree for tree in self.trees])
        consensus_tree = only_trees.get_consensus_tree()
        consensus_tree_with_migrations = TreeWithMigrations(tree=consensus_tree)

        tot_num_trees = len(self.trees)

        migration_events_counter = Counter()
        migration_events_not_in_consensus_tree_counter = Counter()
        num_migration_events_per_tree = None
        for tree_with_migrations in self.trees:
            tree = tree_with_migrations.tree
            migration_events = tree_with_migrations.migration_events
            if num_migration_events_per_tree is None:
                num_migration_events_per_tree = len(migration_events)
            else:
                assert num_migration_events_per_tree == len(migration_events)
            for migration_event in tree_with_migrations.migration_events:
                origin_migration_event_tips = tuple(sorted(tree.get_tip_labels(idx=migration_event.origin_node_id)))
                origin_migration_event_node = consensus_tree.get_mrca_idx_from_tip_labels(origin_migration_event_tips)
                destination_migration_event_tips = tuple(sorted(tree.get_tip_labels(idx=migration_event.destination_node_id)))
                destination_migration_event_node = consensus_tree.get_mrca_idx_from_tip_labels(destination_migration_event_tips)

                if origin_migration_event_node is not None:
                    migration_events_counter[origin_migration_event_node, destination_migration_event_node] += 1
                else:
                    migration_events_not_in_consensus_tree_counter[origin_migration_event_tips, destination_migration_event_tips] += 1

        tot_num_migration_events = sum(migration_events_counter.values()) + sum(migration_events_not_in_consensus_tree_counter.values())

        all_migration_events = []
        for (origin_node_id, destination_node_id), counts in migration_events_counter.items():
            support = round(counts * num_migration_events_per_tree * 100 / tot_num_migration_events)
            migration_event = MigrationEvent(origin_node_id=origin_node_id,
                                             destination_node_id=destination_node_id,
                                             support=support)
            all_migration_events.append(migration_event)

        all_migration_events = list(reversed(sorted(all_migration_events, key=lambda x: x.support)))

        migration_events_not_in_consensus_tree_freqs = {migration_event_ids: round(counts * num_migration_events_per_tree * 100 / tot_num_migration_events) for migration_event_ids, counts in migration_events_not_in_consensus_tree_counter.items()}

        consensus_tree_with_migrations = TreeWithMigrations(consensus_tree,
                                                            migration_events=all_migration_events[:num_migration_events_per_tree])

        return {'consensus_tree': consensus_tree_with_migrations,
                'migration_events_not_in_consensus_tree_freqs': migration_events_not_in_consensus_tree_freqs,
                'all_migration_events': all_migration_events}


def parse_newick_tree(newick_str):
    return toytree.tree(newick_str)


def get_node_idx_from_node_str(node_str, tree):
    if node_str.count(':') == 1:
        tip_labels = [node_str.split(':')[0]]
    else:
        tip_labels = parse_newick_tree(node_str + ';').get_tip_labels()
    return tree.get_mrca_idx_from_tip_labels(tip_labels)


def parse_treemix_tree(treemix):
    lines = treemix.splitlines()
    tree = parse_newick_tree(lines[0])
    
    migration_events = []
    for migration_event_str in lines[1:]:
        items = migration_event_str.split()
        weight = float(items[0])
        origin_node_id = get_node_idx_from_node_str(items[-2], tree)
        destination_node_id = get_node_idx_from_node_str(items[-1], tree)
        migration_event = MigrationEvent(weight=weight,
                                         origin_node_id=origin_node_id,
                                         destination_node_id=destination_node_id)
        migration_events.append(migration_event)
    return TreeWithMigrations(tree=tree, migration_events=migration_events)


def draw_tree(tree, axes, tip_label_mapper=None, tip_label_colors=None,
              tip_label_style=None, draw_support=False, tip_labels=True,
              layout='r', edge_colors=None):

    default_tip_label_style = {#"fill": "#262626",
                               "font-size": "20px",
                               #"-toyplot-anchor-shift": "15px",
                              }
    default_node_labels_style = {"font-size": "20px"}
    default_node_markers = "r2x1.25"
    default_node_sizes = 15

    if tip_label_style is None:
        tip_label_style = default_tip_label_style

    default_tip_color = "#262626"

    if tip_label_mapper is None:
        tip_label_mapper = {}

    if tip_label_colors is None:
        tip_label_colors = {}

    if tip_labels:
        tip_labels = [tip_label_mapper.get(tip_label, tip_label) for tip_label in tree.get_tip_labels()]

    tip_colors = [tip_label_colors.get(tip_label, default_tip_color) for tip_label in tree.get_tip_labels()]

    kwargs = {'axes': axes,
              'tip_labels': tip_labels,
              'tip_labels_colors': tip_colors,
              'tip_labels_style': tip_label_style,
              'use_edge_lengths': True,
              'layout': layout
              }

    if edge_colors:
        kwargs['edge_colors'] =  tree.get_edge_values_from_dict(edge_colors)

    if draw_support:
        kwargs['node_labels'] = 'support'
        #kwargs['node_labels_style'] = default_node_labels_style,
        kwargs['node_markers'] = default_node_markers
        kwargs['node_sizes'] = default_node_sizes

    tree.draw(**kwargs)


def draw_tree_with_migrations(tree_with_migrations, axes,
                              migration_edge_width_multiplier=8,
                              migration_event_colors=None,
                              draw_support=False):

    default_node_labels_style = {"font-size": "20px"}
    default_node_markers = "r2x1.25"
    default_node_sizes = 15

    default_migration_event_color = '#440055'

    #tree_with_migrations.tree._coords.update()

    vert_coords = tree_with_migrations.tree._coords.coords
    vert_coords = tree_with_migrations.tree._coords.verts
    #vert_coords = tree_with_migrations.tree._coords.lines

    draw_tree(tree_with_migrations.tree, axes, draw_support=draw_support)

    if not tree_with_migrations.migration_events:
        return
    
    weights = [migration_event.weight for migration_event in tree_with_migrations.migration_events]
    no_weights = any(weight is None for weight in weights)
    
    if not no_weights:
        min_weight = min(weights)
        max_weight = max(weights)

    if migration_event_colors is None:
        migration_event_colors = [default_migration_event_color] * len(tree_with_migrations.migration_events)

    bootstrap_markers = {'markers': [], 'x_values': [], 'y_values': []}
    for idx, migration_event in enumerate(tree_with_migrations.migration_events):
        origin_coord = vert_coords[migration_event.origin_node_id]
        dest_coord = vert_coords[migration_event.destination_node_id]

        if no_weights:
            width = migration_edge_width_multiplier * 0.5
        else:
            width = migration_event.weight * migration_edge_width_multiplier

        axes.graph(
                numpy.array([[0, 1]]),
                vcoordinates=numpy.array([[origin_coord[0], origin_coord[1]],
                                          [dest_coord[0], dest_coord[1]]]),
                tmarker=">",
                ewidth=width,
                eopacity=0.8,
                vlshow=False,
                vsize=0,
                ecolor=migration_event_colors[idx]
            )   

        support = migration_event.support
        if support is not None:
            bootstrap_markers['x_values'].append((origin_coord[0] + dest_coord[0]) / 2)
            bootstrap_markers['y_values'].append((origin_coord[1] + dest_coord[1]) / 2)
            marker = toyplot.marker.create(shape=default_node_markers, 
                                           label=str(support),
                                           size=default_node_sizes,
                                           mstyle=default_node_labels_style,
                                           #lstyle=nlstyle,
                                          )
            bootstrap_markers['markers'].append(marker)
    if bootstrap_markers:
        axes.scatterplot(bootstrap_markers['x_values'], bootstrap_markers['y_values'],
                         marker=bootstrap_markers['markers'])


def test_parse_newick_tree():
    newick = '(sp_pe:0.0672594,((slc_ec:0.0504711,sp_ec:0.0118801):0.0186141,(slc_ma:0.00382106,(slc_pe:0.00904427,sll_mx:0.0170229):0.00822666):0.152122):0.0443781);'
    tree = parse_newick_tree(newick)
    expected = ['slc_ec', 'slc_ma', 'slc_pe', 'sll_mx', 'sp_ec', 'sp_pe']
    assert sorted(tree.get_tip_labels()) == expected


def test_parse_treemix_tree():
    treemix = '''(sp_pe:0.0672594,((slc_ec:0.0504711,sp_ec:0.0118801):0.0186141,(slc_ma:0.00382106,(slc_pe:0.00904427,sll_mx:0.0170229):0.00822666):0.152122):0.0443781);
0.203407 NA NA NA sp_pe:0.0672594 slc_pe:0.00904427
0.0890222 NA NA NA sp_pe:0.0672594 slc_ma:0.00382106
0.475387 NA NA NA (slc_pe:0.00904427,sll_mx:0.0170229):0.00822666 slc_ec:0.0504711'''
    tree = parse_treemix_tree(treemix)
    assert len(tree.migration_events) == 3
    expected = ['slc_ec', 'slc_ma', 'slc_pe', 'sll_mx', 'sp_ec', 'sp_pe']
    assert sorted(tree.tree.get_tip_labels()) == expected


def test_draw_tree():
    treemix = '''(sp_pe:0.0672594,((slc_ec:0.0504711,sp_ec:0.0118801):0.0186141,(slc_ma:0.00382106,(slc_pe:0.00904427,sll_mx:0.0170229):0.00822666):0.152122):0.0443781);
0.203407 NA NA NA sp_pe:0.0672594 slc_pe:0.00904427
0.0890222 NA NA NA sp_pe:0.0672594 slc_ma:0.00382106
0.475387 NA NA NA (slc_pe:0.00904427,sll_mx:0.0170229):0.00822666 slc_ec:0.0504711'''
    tree_with_migrations = parse_treemix_tree(treemix)
    
    canvas = toyplot.Canvas(style={"background-color": "white"})
    axes = canvas.cartesian()

    #draw_tree(tree_with_migrations.tree, axes)
    draw_tree_with_migrations(tree_with_migrations, axes)
    toyplot.svg.render(canvas, "tree-plot.svg")
    toyplot.png.render(canvas, "tree-plot.png")


def treemix_list():
    treemix1 = '''(sp_pe:0.0672594,((slc_ec:0.0504711,sp_ec:0.0118801):0.0186141,(slc_ma:0.00382106,(slc_pe:0.00904427,sll_mx:0.0170229):0.00822666):0.152122):0.0443781);
0.203407 NA NA NA sp_pe:0.0672594 slc_pe:0.00904427
0.0890222 NA NA NA sp_pe:0.0672594 slc_ma:0.00382106
0.475387 NA NA NA (slc_pe:0.00904427,sll_mx:0.0170229):0.00822666 slc_ec:0.0504711'''
    tree_with_migrations1 = parse_treemix_tree(treemix1)

    treemix2 = '''(sp_pe:0.0672594,((slc_ec:0.0504711,sp_ec:0.0118801):0.0186141,(slc_ma:0.00382106,(slc_pe:0.00904427,sll_mx:0.0170229):0.00822666):0.152122):0.0443781);
0.203407 NA NA NA sp_pe:0.0672594 slc_pe:0.00904427
0.0890222 NA NA NA sp_pe:0.0672594 sll_mx:0.00382106
0.475387 NA NA NA (slc_pe:0.00904427,sll_mx:0.0170229):0.00822666 slc_ec:0.0504711'''
    tree_with_migrations2 = parse_treemix_tree(treemix2)


    trees = TreeWithMigrationsList([tree_with_migrations1] * 100 + [tree_with_migrations2] * 100)

    canvas = toyplot.Canvas(style={"background-color": "white"})
    axes = canvas.cartesian()
    consensus_tree = trees.get_consensus_tree()['consensus_tree']
    consensus_tree.tree = consensus_tree.tree.root(['sp_pe'])
    draw_tree_with_migrations(consensus_tree,
                              axes,
                              draw_support=True)
    axes.show = False
    toyplot.png.render(canvas, "tree-plot-consensus.png")

if __name__ == '__main__':
    test_parse_newick_tree()
    test_parse_treemix_tree()
    test_draw_tree()
    treemix_list()
    