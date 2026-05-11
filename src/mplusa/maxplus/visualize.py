import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from .geometry import AbstractPolytope, project_point


def hasse_diagram(polytope : AbstractPolytope) -> None:
    """ Draws and shows the Hasse diagram of the given polytope using matplotlib. """
    _, ax = plt.subplots()
    ax.set_xticks([])
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')
    ax.set_ylabel('Rank')
    for spine in ['top', 'left', 'bottom']:
        ax.spines[spine].set_visible(False)
    width = max(map(len, polytope.structure.values()))
    layers = {}
    for rank, layer in polytope.structure.items():
        points = []
        for i, connections in enumerate(layer):
            x = (i * width) / (len(layer) - 1) if len(layer) > 1 else width / 2
            points.append((x, rank))
            if rank > 0:
                for connection in connections:
                    ax.plot(
                        [x, layers[rank - 1][connection][0]],
                        [rank, layers[rank - 1][connection][1]],
                        color='black'
                    )
        layers[rank] = points
        for i, point in enumerate(points):
            ax.text(
                point[0], point[1], str(i),
                ha='center', va='center', fontsize=10,
                bbox={
                    'facecolor': 'white',
                    'edgecolor': 'black',
                    'boxstyle': 'circle',
                }
            )
    plt.show()


def draw_polytope2D(
    polytope : AbstractPolytope,
    color : str = 'black',
    show_vertices : bool = True,
    vertices_color : str|None = None,
    vertices_marker : str = 'o',
    show_pseudovertices : bool = False,
    pseudovertices_color : str|None = None,
    pseudovertices_marker : str = 'o',
    show : bool = False) -> None:
    """ Draws a given polytope on a 2-dimensional surface using matplotlib. """
    if polytope.dimension < 2:
        raise NotImplementedError('Projection of polytopes of lesser dimensions not implemented currently.')
    if polytope.dimension != 2:
        raise ValueError('The polytope cannot be projected onto a 2-dimensional surface.')
    if show_vertices:
        vertices = list(map(project_point, polytope.vertices))
        if vertices_color is None:
            vertices_color = color
        plt.scatter(
            [vertex[0] for vertex in vertices],
            [vertex[1] for vertex in vertices],
            color=vertices_color,
            marker=vertices_marker
        )
    if show_pseudovertices:
        pseudovertices = list(map(project_point, polytope.pseudovertices))
        if pseudovertices_color is None:
            pseudovertices_color = color
        plt.scatter(
            [pseudovertex[0] for pseudovertex in pseudovertices],
            [pseudovertex[1] for pseudovertex in pseudovertices],
            color=pseudovertices_color,
            marker=pseudovertices_marker
        )
    line_segments = polytope.get_all_line_segments()
    for line_segment in line_segments:
        line_segment = list(map(project_point, line_segment))
        plt.plot(
            [vertex[0] for vertex in line_segment],
            [vertex[1] for vertex in line_segment],
            color=color
        )
    if show:
        plt.show()


def draw_hyperplane2D(
    hyperplane : Hyperplane,
    color : str = 'black',
    scale : float = 1,
    show : bool = False) -> None:
    """ Draws a given hyperplane on a 2-dimensional surface using matplotlib. """
    apex = project_point(hyperplane.get_apex())
    plt.plot([apex[0] - scale, apex[0]], [apex[1], apex[1]], color=color)
    plt.plot([apex[0] + scale, apex[0]], [apex[1] + scale, apex[1]], color=color)
    plt.plot([apex[0], apex[0]], [apex[1] - scale, apex[1]], color=color)
    if show:
        plt.show()
