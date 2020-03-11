
def get_sorted_legend_handles(axes):

    handles, labels = axes.get_legend_handles_labels()
    # sort both labels and handles by labels
    handles, labels = zip(*sorted(zip(handles, labels), key=lambda t: t[1]))
    return handles, labels
