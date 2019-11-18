import argparse

from bokeh.layouts import row, column, widgetbox
from bokeh.models import HoverTool, Legend, LegendItem, Select
from bokeh.plotting import figure, curdoc, ColumnDataSource

from woldrnaseq.models import (
    load_all_coverage,
    load_experiments,
    load_library_tables,
)


class Coverage:
    """Common code for coverage plots
    """
    def __init__(self, experiments, libraries):
        self.experiments = experiments
        self.libraries = libraries
        self.experiment_names = sorted(self.experiments.index)

        self.coverage = load_all_coverage(self.libraries)

        self.experiments_combo = Select(
            title='Experiments',
            value=self.experiment_names[0],
            options=self.experiment_names)
        self.experiments_combo.on_change('value', self.update_experiment)

        self._layout = None

    @property
    def experiment_name(self):
        """Return name of the currently selected experiment
        """
        return self.experiments_combo.value

    @experiment_name.setter
    def experiment_name(self, value):
        self.experiments_combo.value = value

    @property
    def library_names(self):
        """Return list of library IDs for the currently selected experiment
        """
        return self.experiments.loc[self.experiment_name].replicates

    def use_experiment(self, value):
        if value is None:
            pass
        elif value in self.experiment_names:
            self.experiment_name = value
        else:
            raise ValueError("{} not a valid experiment name".format(value))

    def make_plot(self, experiment_name, library_id=None):
        raise NotImplementedError("This is a base class")

    def static_layout(self):
        self._layout = self.make_plot(self.experiment_name)
        return self._layout
