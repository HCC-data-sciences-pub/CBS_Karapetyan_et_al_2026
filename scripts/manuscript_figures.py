import itertools
from pathlib import Path
from typing import *
from hcc_164_scripts import projectinfo
import pandas
from pprint import pprint
import matplotlib.pyplot as plt
import seaborn
from loguru import logger
from scipy import stats
import math
from vectratools import utilities

plt.style.use('fivethirtyeight')

LEADIN_ORDER = ['Rela', 'Nivo', 'Combo']


def get_compartment_tissue_type(id_compartment: str) -> Literal['normal', 'tumor']:
	suffix = id_compartment.split('.')[0].split('_')[-1]
	result = 'normal' if suffix.startswith('N') else 'tumor'
	return result


def make_boxplot(table_pairplot: pandas.DataFrame, marker_column: str = None, path_output: Path = None, title: str = None):
	palette_leadin = {
		'Combo': '#ee3bee',
		'Nivo':  '#00cdcd',
		'Rela':  '#737373'
	}
	color_background = "#FFFFFF"

	options = {
		# 'pvalueType':       "comparison",  # Literal['clustering', 'comparison'],
		'pvalueType':       'comparison',
		'color:background': '#FFFFFF',
		'palette:leadin':   {},
		'palette:dataset':  {},
	}

	_base_figsize = 10
	aspect_ratio = 0.73  # height-to-width

	fontsize_title = 24
	fontsize_axis_label_x = 24
	fontsize_xticklabels = 18
	fontsize_axis_label_y = 18
	fontsize_text = 18

	table_pairplot_baseline = table_pairplot

	y_min = table_pairplot_baseline[marker_column].min()
	y_max = table_pairplot_baseline[marker_column].max()
	y_max_abs = max(abs(y_min), abs(y_max)) + 1000  # Add a small buffer so the maximal point isn't cut off.
	# fig, ax = plt.subplots(figsize = (10, 15))
	fig, ax = plt.subplots(figsize = (_base_figsize, _base_figsize * aspect_ratio))

	ax.axhline(0, color = 'k', linestyle = 'dotted')

	ax = seaborn.boxplot(
		table_pairplot_baseline,
		x = 'leadin', y = marker_column, hue = 'leadin',
		order = LEADIN_ORDER,
		palette = palette_leadin,
		ax = ax,
		fill = False,
		fliersize = 0
	)

	ax = seaborn.stripplot(
		table_pairplot_baseline,
		x = 'leadin', y = marker_column, hue = 'leadin',
		order = LEADIN_ORDER,
		palette = palette_leadin,
		edgecolor = 'k',
		linewidth = 2,
		size = 10,
		ax = ax,
		dodge = False
	)

	ax.set_ylim(-y_max_abs * 1.1, y_max_abs * 1.1)
	ax.set_ylim(-21000, 110000)

	if ax.get_legend() is not None:
		ax.get_legend().remove()

	#################################### Set Plot Labels #######################################
	if title is None:
		title = "All TMA and Tumor Wholeslide Samples\nBaseline Samples"
	fig.suptitle(title, fontsize = fontsize_title)

	ax.set_ylabel(f"Ripley's K\n{marker_column.replace('_', ' vs ')}", fontsize = fontsize_axis_label_y)
	ax.set_xlabel("Lead-in arms", fontsize = fontsize_axis_label_x)
	ax.tick_params(axis = 'x', which = 'major', labelsize = 24)

	############ Change the color of the plot background ###################
	# Set background color of the figure
	fig.patch.set_facecolor(color_background)
	plt.rcParams['figure.facecolor'] = color_background

	# Set background color of the subplot
	ax.set_facecolor(color_background)

	################ Change color of the plot frame ###########################
	ax.spines['bottom'].set_color(color_background)
	ax.spines['top'].set_color(color_background)
	ax.spines['left'].set_color(color_background)
	ax.spines['right'].set_color(color_background)

	# ax.text(-1 * (10/_base_figsize), y_max_abs / 2, "Clustering", fontsize = fontsize_text, rotation = 90)
	# ax.text(-1*(10/_base_figsize), -y_max_abs / 2 - (y_max_abs / 4), "Dispersion", fontsize = fontsize_text, rotation = 90)
	ax.text(-1 * (10 / _base_figsize) + 0.10, 10000, "Clustering", fontsize = fontsize_text, rotation = 90, color = 'tab:cyan')
	ax.text(-1 * (10 / _base_figsize) + 0.10, -30000, "Dispersion", fontsize = fontsize_text, rotation = 90, color = 'tab:red')

	ax = add_pvalues(
		table_pairplot_baseline,
		ax,
		pvalue_type = options['pvalueType'],
		marker_column = marker_column
	)

	plt.tight_layout()
	if path_output:
		plt.savefig(path_output.with_suffix('.png'), facecolor = color_background)
		plt.savefig(path_output.with_suffix('.svg'), facecolor = color_background)
	plt.show()


def add_pvalues(table_pairplot_baseline: pandas.DataFrame, ax: plt.Axes, pvalue_type: Literal['clustering', 'comparison'], marker_column: str):
	series_combo = table_pairplot_baseline[table_pairplot_baseline['leadin'] == 'Combo'][marker_column]
	series_nivo = table_pairplot_baseline[table_pairplot_baseline['leadin'] == 'Nivo'][marker_column]
	series_rela = table_pairplot_baseline[table_pairplot_baseline['leadin'] == 'Rela'][marker_column]

	median_is_positive = series_nivo.median() > 0
	y_value_max = max([series_nivo.max(), series_rela.max(), series_combo.max()])
	y_value_min = min([series_nivo.min(), series_rela.min(), series_combo.min()])
	y_value_largest = max(abs(y_value_min), abs(y_value_max))

	if median_is_positive:
		# y_values = [series_nivo.min(), series_rela.min(), series_combo.min()]
		y_value = (y_value_min - 5000) * 1.1
	else:
		# y_values = [series_nivo.max(), series_rela.max(), series_combo.max()]
		y_value = y_value_max * 1.1

	if abs(y_value) < (y_value_largest / 2):
		y_value = math.copysign(y_value_largest / 2, y_value)
	elif abs(y_value) * 1.2 > y_value_largest:
		y_value = math.copysign(y_value_largest * 0.80, y_value)

	y_value = 40000

	if pvalue_type == 'clustering':
		alternative = 'greater' if median_is_positive else 'less'
		p_value_nivo = stats.quantile_test(series_nivo.values, alternative = alternative)
		p_value_rela = stats.quantile_test(series_rela.values, alternative = alternative)
		p_value_combo = stats.quantile_test(series_combo.values, alternative = alternative)

		x_offset = 0.25
		fontsize = 24
		ax.text(0 - x_offset, y_value, f"p={p_value_nivo.pvalue:.3f}", fontsize = fontsize)
		ax.text(1 - x_offset, y_value, f"p={p_value_rela.pvalue:.3f}", fontsize = fontsize)
		ax.text(2 - x_offset, y_value, f"p={p_value_combo.pvalue:.3f}", fontsize = fontsize)
	else:
		use_ranksums = True
		if use_ranksums:
			p_value_nivo_rela = stats.ranksums(series_nivo.values, series_rela.values).pvalue
			p_value_nivo_combo = stats.ranksums(series_nivo.values, series_combo.values).pvalue
			p_value_rela_combo = stats.ranksums(series_rela.values, series_combo.values).pvalue
		else:
			result = stats.tukey_hsd(series_nivo.values, series_rela.values, series_combo.values)
			p_value_nivo_rela = result.pvalue[0][1]
			p_value_nivo_combo = result.pvalue[0][2]
			p_value_rela_combo = result.pvalue[1][2]

		add_whisker(ax, (0, 1), y_value, f"{p_value_nivo_rela:.03f}")
		# add_whisker(ax, (0, 2), y_value * 1.1, f"{p_value_nivo_combo:.03f}")
		add_whisker(ax, (0, 2), y_value * 1.2, f"{p_value_rela_combo:.03f}")
		# add_whisker(ax, (1, 2), y_value * 1.2, f"{p_value_rela_combo:.03f}")
		add_whisker(ax, (1, 2), y_value * 1.4, f"{p_value_nivo_combo:.03f}")

	return ax


def add_whisker(ax: plt.Axes, xlims: Tuple[float, float], y_value: float, label: str, ratio: float = 0.02) -> plt.Axes:
	color_whisker = 'k'
	linewidth = 1
	whiskery = (y_value - (y_value * ratio), y_value + (y_value * ratio))

	# Add the horizontal part
	ax.plot(xlims, [y_value, y_value], color = color_whisker, linewidth = linewidth)

	# Add the vertical parts
	ax.plot([xlims[0], xlims[0]], whiskery, color = color_whisker, linewidth = linewidth)
	ax.plot([xlims[1], xlims[1]], whiskery, color = color_whisker, linewidth = linewidth)

	# Add the label
	y_offset = y_value * 0.05
	if y_offset < 0:
		y_offset *= 2
	ax.text((sum(xlims) / 2) - 0.1, y_value + y_offset, label)

	return ax


def merge_compartments(table_pairplot: pandas.DataFrame) -> pandas.DataFrame:
	"""
		Merges TMA and Vectra samples from the same sample.
	"""
	columns = [i for i in table_pairplot.columns if '-' in i]
	table_pairplot.loc[:, "id:sample"] = table_pairplot['id:compartment'].apply(lambda s: s.rpartition("_")[0])
	if 'tissueType' not in table_pairplot.columns:
		table_pairplot.loc[:, 'tissueType'] = [get_compartment_tissue_type(i) for i in table_pairplot['id:compartment'].tolist()]

	groups = table_pairplot.groupby(by = ['id:sample', 'id:dataset', 'tissueType'])

	rows = list()
	for (id_sample, id_dataset, id_tissuetype), group in groups:
		row = group.iloc[0]
		if len(group) == 1:
			rows.append(row.to_dict())
		else:
			record = {key: row[key] for key in ['id:sample', 'id:dataset', 'tissueType', 'leadin', 'response:P04', 'response:P16', 'response:best', 'timepoint']}
			group_mean = group[columns].mean()
			record |= group_mean.to_dict()
			rows.append(record)
	logger.debug(f"{len(table_pairplot)=}")
	df = pandas.DataFrame(rows)
	logger.debug(f"{len(df)=}")

	return df


def main():
	folder_project = Path.home() / "Documents" / "projects" / "HCC-CBS-164-Hillman-TBruno"

	options = {
		'outputFormat':  'png',  # Literal['png', 'pdf']
		'title':         None,
		'tissueType':    'normal',  # Literal['normal', 'tumor']
		'timepoint':     "baseline",  # Literal['baseline', 'week 4']'
		'datasetLabel':  'source',
		'input:table':   folder_project / "workflows" / "workflow_source_2024-02-29" / "ripley_bivariate_values.source.tsv",
		'output:folder': Path("/home/proginoskes/Documents/projects/HCC-CBS-164-Hillman-TBruno/reports/2026-10-16_report/figures")
		# 'input:table': folder_project / "workflows" / "workflow_LAG3_2025-03-15_tiles" /"spatialtime" / "ripley_bivariate_values.50um.tsv",
		# 'output:folder': utilities.checkdir(folder_project / "workflows" / "workflow_LAG3_2025-03-15_tiles" / "spatialtime" / "figures" / "boxplots")
	}
	path_output = folder_project / "reports" / "2026-10-16_report" / "CD8_FOXP3.ripley.BL.Tumor.png"

	############################################################ Read in the input table ###################################################################
	table_pairplot = pandas.read_csv(options['input:table'], sep = "\t")

	# Add a few new columns to help with filtering the dataset.
	table_pairplot.loc[:, 'which:dataset'] = options['datasetLabel']
	table_pairplot.loc[:, 'id:patient'] = table_pairplot['id:compartment'].apply(lambda s: s.split('_')[0])
	# The tissue type is determind by the compartment barcode.
	table_pairplot.loc[:, 'tissueType'] = [get_compartment_tissue_type(i) for i in table_pairplot['id:compartment'].tolist()]
	table_pairplot.columns = [i.replace(':', '_').replace('-', '_') for i in table_pairplot.columns]

	# 740 x 540
	title = f"CD8_vs_FOXP3 - Tumor - BL"
	table_pairplot_timepoint = table_pairplot[(table_pairplot['tissueType'] == 'tumor') & (table_pairplot['timepoint'] == 'BL')]
	make_boxplot(table_pairplot_timepoint, marker_column = 'CD8_FOXP3', path_output = path_output, title = title)


def main_aggregate_ripleys_k():
	folder_tables = Path("/home/proginoskes/archive/workspace/HCC-CBS-164-Hillman-TBruno/workflows/workflow_source_2024-02-29/spatialtime/results_exact")
	folder_compartments = Path("/home/proginoskes/archive/workspace/HCC-CBS-164-Hillman-TBruno/workflows/workflow_source_2024-02-29/spatialtime/compartments")
	paths_compartments = {path.stem.split('.')[0]: path for path in folder_compartments.glob("*.tsv")}

	results = list()
	for path in folder_tables.glob("*bivariate.tsv"):

		df_compartment = pandas.read_csv(paths_compartments[path.stem.split('.')[0]], sep = "\t")
		first_row = df_compartment.iloc[0]

		df = pandas.read_csv(path, sep = "\t")

		groups = df.groupby(by = ['Anchor', 'Counted'])

		for (anchor_column, counted_column), group in groups:
			row = group[group['r'] == 50].iloc[0]

			record = {
				'barcode':            row['id:mask'],
				'timepoint':          first_row['id:timepoint'],
				'response':           first_row['response:overallBest'],
				'anchor':             anchor_column.split(':')[-1],
				'counted':            counted_column.split(':')[-1],
				'degreeOfClustering': row['Degree of Clustering Exact']
			}
			results.append(record)

	df = pandas.DataFrame(results)
	return df

	path = Path("/home/proginoskes/Documents/projects/HCC-CBS-164-Hillman-TBruno/workflows/workflow_source_2024-02-29/ripley_bivariate_values.source.tsv")
	table = pandas.read_csv(path, sep = "\t")
	df = merge_compartments(table)
	print(len(table))
	print(len(df))


if __name__ == "__main__":
	main()
