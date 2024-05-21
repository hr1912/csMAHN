#!/usr/bin/env python
# coding: utf-8

# ```shell
# conda activate
# cd ~/link/res_publish/plot
# jupyter nbconvert init_plot.ipynb --to python
# echo 'finish'
# 
# 
# ```

# In[1]:


import sys
from pathlib import Path
p_root = Path('~/link/res_publish').expanduser()
None if str(p_root) in sys.path else sys.path.append(str(p_root))


# In[2]:


from func import *

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

# cm = mpl.colormaps['magma']
# cm = mpl.colormaps['Reds']
cm = mpl.colormaps['YlOrRd']
cm_2 = mpl.colormaps['magma']


# # plot para
# 

# In[3]:


ppara_adata_umap = {}
ppara_color_map = {}
ppara_data = {}
ppara_plot_function_custom = {}


# In[4]:


mpl.font_manager.fontManager.addfont(p_plot.joinpath('font', 'arial.ttf'))
# print(*mpl.font_manager.get_font_names(),sep='\n')
sns.set_theme(rc=mpl.rcParams)
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['font.family'] = 'arial'
mpl.rcParams['font.size'] = 6  # 磅（points）
mpl.rcParams['axes.facecolor'] = 'white'
mpl.rcParams['legend.title_fontsize'] = mpl.rcParams['font.size']
mpl.rcParams['axes.titlesize'] = mpl.rcParams['font.size'] + 2
mpl.rcParams['axes.titleweight'] = 'bold'
mpl.rcParams['axes.labelsize'] = mpl.rcParams['font.size']
mpl.rcParams['xtick.labelsize'] = mpl.rcParams['font.size']
mpl.rcParams['ytick.labelsize'] = mpl.rcParams['font.size']
mpl.rcParams['patch.edgecolor'] = 'white'
mpl.rcParams['patch.linewidth'] = .5

fontdict_default = {
    _k: mpl.rcParams[_k.replace('font', 'font.')]
    for _k in 'fontfamily,fontsize,fontweight'.split(',')}
fontdict_axes_title = fontdict_default.copy()
fontdict_axes_title.update({'fontsize': mpl.rcParams['axes.titlesize']})


# # plt_pdf
# 
# 使用matplot绘制出A4大小的fig,在保存为pdf,避免了后期使用AI的组图（AI组合多个pdf是在是太卡了）
# 
# + `plt_pdf_add_grid_customer`将fig分为多个网格,每个网格的长宽比接近1:1
# 
# + `plt_pdf_add_grid_customer`显示该网格
# 
# + `plt_pdf_add_ax_with_spec`通过选择不同的网格坐标和长宽决定所绘制的图在fig上的位置和大小
#     + 当需要迁移至其他页仅需更改fig
#     + 当需要移动位置，仅需更改x,y
#     + 当需要ax大小，仅需更改x_offest,y_offest
# 
# + `plt_pdf_add_text_with_ax`添加文字

# In[5]:


def plt_pdf_add_grid_fig_coordinate(fig, alpha=0.35):
    """标出fig的绝对坐标系
"""
    for i, (v) in enumerate(np.linspace(0, 1, 10+1)):
        # |
        fig.add_artist(
            Line2D([v, v], [0, 1], linestyle='--', alpha=alpha, c='#2F4F4F'))
        for y in [0.1, 0.5, 0.9]:
            fig.text(
                v+0.01,
                y,
                'x=.{:.0f}'.format(i) if v == 0 or v == .5 else '.{:.0f}'.format(i),
                fontdict={
                    'fontsize': 10,
                    'color': 'black',
                    'alpha': alpha})
    fontdict = {'fontsize': 10, 'color': 'blue', 'alpha': alpha}

    for i, (v) in enumerate(np.linspace(0, 1, 10 + 1)):
        # ------
        fig.add_artist(
            Line2D([0, 1], [v, v], linestyle='--', alpha=alpha, c='#FFA500'))
        for x in [0.1, 0.5, 0.9]:
            fig.text(
                x,
                v+0.01,
                'y=.{:.0f}'.format(i) if v == 0 or v == .5 else '.{:.0f}'.format(i),
                fontdict={
                    'fontsize': 10,
                    'color': 'orange',
                    'alpha': alpha})


def plt_pdf_add_grid_customer(fig, alpha=0.15):
    """为fig在以下参数下进行了优化
display([8.27*.97, 11.96*.96])
display([8.27*.97/.25, 11.96*.96/.25])
# 以0.25为单位1，将width分为32等分，higth分为46等分
fig = plt.figure(figsize=(8.27, 11.69))
spec = fig.add_gridspec(nrows=46, ncols=32,
                    left=0.03, right=1,  # 设置边距
                    bottom=0.02, top=0.98,  # 设置边距
                    wspace=0, hspace=0)  # 设置子图间距
"""
    for i, (v) in enumerate(np.linspace(0.03, 1, 32+1)):
        # |
        fig.add_artist(Line2D([v, v], [0.02, 0.98],
                       linestyle='--', alpha=alpha, c='#00BFFF'))
        if i >= 32 or i % 2 == 1:
            continue
        for y in [0.1, 0.5, 0.9]:
            fig.text(
                v+0.01,
                y,
                'x={:.0f}'.format(i) if i == 0 or i == 16 else '{:.0f}'.format(i),
                fontdict={
                    'fontsize': 10,
                    'color': 'blue',
                    'alpha': alpha})
    fontdict = {'fontsize': 10, 'color': 'blue', 'alpha': alpha}
    for i, (v) in enumerate(np.linspace(0.98, 0.02, 46 + 1)):
        # ------
        fig.add_artist(Line2D([0.03, 1], [v, v],
                              linestyle='--', alpha=alpha, c='#D02090'))
        if i >= 46 or i % 2 == 1:
            continue
        for x in [0.1, 0.5, 0.9]:
            fig.text(
                x,
                v-0.01,
                'y={:.0f}'.format(i) if i == 0 or i == 24 else '{:.0f}'.format(i),
                fontdict={
                    'fontsize': 10,
                    'color': 'red',
                    'alpha': alpha})


def plt_pdf_subset_spec(spec, x, y, x_offest=1, y_offest=1):
    return spec[y:y+y_offest, x:x+x_offest]


def plt_pdf_add_ax_with_spec(fig, spec, x, y, x_offest=1, y_offest=1):

    return fig.add_subplot(
        plt_pdf_subset_spec(
            spec, x, y, x_offest, y_offest))


def plt_pdf_add_text_with_ax(
    ax,
    text,
    x=0,
    y=1,
    fontdict={
        'fontsize': 14,
        'fontweight': 'bold'}):
    ax.text(x, y, text, fontdict=fontdict)
    ax.set_axis_off()


# # Block
# 
# `块` 用于将代码分块
# 
# 因为在notebook中用注释将代码分块，无法实现代码的折叠，使得代码较为混乱
# 
# 故通过`with Block():`组合构造出可以折叠的with语句块，进而提高代码的可读性
# 
# + `with Block():`内部并未与外部进行隔离
# 
# + 实现了计时功能

# In[6]:


class Block:
    """用于在notebook中给代码划分区域(块),从而使代码能够折叠
with Block('test',show_comment=True):
    print('inner')
    time.sleep(2)
# output
## [start][test] 0509-00:20:47------------------------------------------------
## inner
## [end][test] 0509-00:20:49--------------------------------------------------
        2.002 s used
    """

    def __init__(self, comment, show_comment=False):
        self.comment = comment
        self.show_comment = show_comment
        self.time = 0

    def __enter__(self):
        if self.show_comment:
            self.time = time.time()
            print("[start][{}] {}".format(
                self.comment,
                time.strftime('%m%d-%H:%M:%S',
                              time.localtime())).ljust(75, '-'))

    def __exit__(self, exc_type, exc_value, traceback):
        if self.show_comment:
            print(
                "[end][{}] {}".format(
                    self.comment,
                    time.strftime(
                        '%m%d-%H:%M:%S',
                        time.localtime())).ljust(
                    75,
                    '-'))
            print("\t{:.3f} s used".format(time.time()-self.time))


# # EnrichmentAnalysis
# 
# 借由脚本`Rscript EnrichmentAnalysis.r` (R 富集分析器)实现在python通过shell调用指定conda环境中的R
# 
# 进行EnrichAnalysis

# In[15]:





# In[7]:


def EnrichmentAnalysis_run_with_r(
        tag,
        geneset_label,
        gene_string,
        separator=',',
        env='publish'):
    command = "{}Rscript EnrichmentAnalysis.r 'data/EnrichmentAnalysis/{}' '{}' '{}' {}".format(
        '~/apps/miniconda3/envs/{}/bin/'.format(env) if env else '',
        tag, geneset_label, gene_string, separator

    )
    os.system(command)


def EnrichmentAnalysis_show_genesets(p=Path(
        '~/link/res_publish/plot/data/EnrichmentAnalysis').expanduser()):
    gs_info = pd.read_csv(p.joinpath('geneset/geneset_info.csv'))
    display(gs_info.loc[:, 'label,version'.split(',')])


def EnrichmentAnalysis_find_genesets_path(geneset_label, p=Path(
        '~/link/res_publish/plot/data/EnrichmentAnalysis').expanduser()):
    gs_info = pd.read_csv(p.joinpath('geneset/geneset_info.csv'))
    gs_info['path'] = gs_info['path'].apply(
        lambda x: p.joinpath('geneset', x))
    res = gs_info.query("label == '{}'".format(geneset_label))
    if res.shape[0] == 1:
        res = res['path'].to_list()[0]
    else:
        raise Exception("[Error] find {} itme with geneset_label='{}'"
                        .format(res.shape[0], geneset_label))
    return res


def EnrichmentAnalysis_get_res_df(
        p=p_plot.joinpath('data/EnrichmentAnalysis')):
    assert p.exists(), '[not exists] {}'.format(p)
    df = pd.DataFrame({'dir': p.iterdir()})
    df = df[~df['dir'].apply(lambda x: x.match('*/.ipynb_checkpoints'))]
    df.index = df['dir'].apply(lambda x: x.name)
    df['p_gene'] = df['dir'].apply(lambda x: x.joinpath('gene_input.txt'))
    df['p_table'] = df['dir'].apply(
        lambda x: x.joinpath('Enrich_table.csv'))
    df['p_table_all'] = df['dir'].apply(
        lambda x: x.joinpath('Enrich_table_all.csv'))
    return df


def EnrichmentAnalysis_get_adj(row, all=False):
    return pd.read_csv(
        row['p_table{}'.format('_all' if all else '')],
        index_col=0).loc[:, ['p.adjust']].rename(columns={
            'p.adjust': row.name
        })


# # sc_pl_show_genes
# 
# 借由scanpy 进行绘图, 并做了些许调整
# 
# 但是可恶的错误提示仍然未能避免，气

# In[8]:


funcs_select_axes = {
    '0_colorbar': lambda fig,
    marker: [
        i for i in fig.get_axes() if i.get_title() == 'z-score of\nmean expression'],
    '1_box': lambda fig,
    marker: [
        i for i in fig.get_axes() if list(
            marker.keys()) == [
            i1.get_text() for i1 in i.get_yticklabels()]],
}


def sc_pl_stacked_violin(adata, marker, group_by, categories_order=None,
                         ax=None, del_yticks=False, var_group_rotation=0,
                         funcs_select_axes=funcs_select_axes):
    if categories_order is None:
        categories_order = np.sort(adata.obs[group_by].unique())
    else:
        categories_order = [_ for _ in categories_order
                            if _ in adata.obs[group_by].unique()]

    sc.pl.stacked_violin(
        adata, marker, group_by,
        var_group_rotation=var_group_rotation,
        colorbar_title="z-score of\nmean expression",
        vmax=2.5, vmin=-2.5, cmap=cm,  # cmap="Blues",
        # layer="scaled",
        dendrogram=False, swap_axes=False,
        ax=ax, show=False,  # return_fig=True
    )
    fig = ax.figure

    for i in funcs_select_axes['0_colorbar'](fig, marker):
        fig.delaxes(i)
    for i in funcs_select_axes['1_box'](fig, marker):
        i.tick_params(length=0, width=0)
        i.set_frame_on(False)
        i.set_yticks([], []) if del_yticks else None

    for _patches in np.concatenate([i.patches for i in fig.get_axes()
                                    if list(marker.keys()) == [
        i2.get_text() for i2 in i.get_children()
            if isinstance(i2, mpl.text.Text) and len(i2.get_text()) > 0]
    ]):
        _vertices = _patches.get_path().vertices
        _vertices[:, 1] = np.where(
            _vertices[:, 1] == 0.6, 0.1, _vertices[:, 1])
        _patches.set_path(
            mpl.path.Path(
                _vertices,
                _patches.get_path().codes))
        _patches.set(linewidth=.75, linestyle='-', color='grey')
    fig.draw_without_rendering()


def sc_pl_matrixplot(adata, marker, group_by, categories_order=None,
                     ax=None, del_yticks=False, var_group_rotation=0,
                     funcs_select_axes=funcs_select_axes):
    sc.pl.matrixplot(
        adata, marker, group_by,
        var_group_rotation=var_group_rotation,
        colorbar_title="z-score of\nmean expression",
        vmax=2.5, vmin=-2.5, cmap=cm,  # cmap="Blues",
        layer="scaled",
        dendrogram=False, swap_axes=False,
        ax=ax, show=False,  # return_fig=True
    )
    fig = ax.figure

    for i in funcs_select_axes['0_colorbar'](fig, marker):
        fig.delaxes(i)
    for i in funcs_select_axes['1_box'](fig, marker):
        i.tick_params(length=0, width=0)
        i.set_frame_on(False)
        i.set_yticks([], []) if del_yticks else None

    for _patches in np.concatenate([i.patches for i in fig.get_axes()
                                    if list(marker.keys()) == [
        i2.get_text() for i2 in i.get_children()
            if isinstance(i2, mpl.text.Text) and len(i2.get_text()) > 0]
    ]):
        _vertices = _patches.get_path().vertices
        _vertices[:, 1] = np.where(
            _vertices[:, 1] == 0.6, 0.1, _vertices[:, 1])
        _patches.set_path(
            mpl.path.Path(
                _vertices,
                _patches.get_path().codes))
        _patches.set(linewidth=2, linestyle='-', color='grey')
    fig.draw_without_rendering()


map_func_sc_pl_show_genes = {
    'stacked_violin': sc_pl_stacked_violin,
    'matrixplot': sc_pl_matrixplot,
}


def sc_pl_show_genes(key, adata, marker, group_by, categories_order=None,
                     ax=None, del_yticks=False, var_group_rotation=0,
                     funcs_select_axes=funcs_select_axes):
    """
    key:
        stacked_violin
        matrixplot
    """

    map_func_sc_pl_show_genes[key](
        adata, marker=marker, group_by=group_by,
        categories_order=categories_order,
        var_group_rotation=var_group_rotation,
        ax=ax, del_yticks=del_yticks,
        funcs_select_axes=funcs_select_axes
    )


# In[9]:


ppara_data['key_scpl_show_genes'] = 'stacked_violin'


# # heatmap_gene
# 
# 借助sns实现的,仿照sc.pl.matrixplot

# In[ ]:


# with Block("ppara_plot_function_custom['heatmap_gene']"):
#     def _func(
#             df_plot,
#             df_marker,
#             ax,
#             zscore=True,
#             cmap=cm,
#             cbar=False,
#             line_h=.25,
#             kv_line={
#                 'linestyle': '-',
#                 'linewidth': 1,
#                 'color': 'grey',
#                 'dash_joinstyle': 'miter',
#                 'dash_capstyle': 'butt'},
#             kv_heatmap={
#                 'vmax': 2.5,
#                 'vmin': -2.5},
#             **kvargs):
#         """
#            kvargs:
#                del_yticks: default False
#                kv_line_text: default fontdict_default

# """
#         df_marker.index = np.arange(df_marker.shape[0])
#         if df_plot.columns[-1].startswith('gap_'):
#             df_plot = df_plot.iloc[:, :-1].copy()
#         if zscore:
#             df_plot = stats.zscore(df_plot, axis=0,nan_policy='omit')
#         sns.heatmap(df_plot, square=True, cbar=cbar,
#                     cmap=cmap, ax=ax, **kv_heatmap)
#         if kvargs.setdefault('del_yticks', False):
#             ax.set_yticks([], [], rotation=0,
#                           **fontdict_default)
#         else:
#             _texts = df_plot.index
#             ax.set_yticks(
#                 np.arange(
#                     df_plot.shape[0])+.5,
#                 df_plot.index,
#                 rotation=0,
#                 **fontdict_default)
#         _texts = df_plot.columns.str.replace('^gap_\\d+$', '', regex=True)
#         ax.set_xticks(np.arange(df_plot.shape[1])+.5, _texts, rotation=90,
#                       **fontdict_default)
#         # _temp = [-1,0,1]
#         # ax.set_yticks(_temp,[str(_) for _ in _temp])
#         ax.set_ylim(min(ax.get_ylim()), max(ax.get_ylim())+1)

#         # plot line
#         for text, s, e in zip(
#                 df_marker.query('line_start')['cell_type'],
#                 df_marker.query('line_start').index,
#                 df_marker.query('line_end').index.to_list() + [df_marker.index.size]):
#             h_s, h = df_plot.shape[0] + .25, line_h
#             s = s+.25
#             e = e-.25
#             # text = 'mast cell'
#             ax.step([s, s, e, e], [h_s, h_s+h, h_s+h, h_s], **kv_line)
#             _fontdict = fontdict_default.copy()
#             _fontdict.update({'horizontalalignment': 'center',
#                               'multialignment': 'center',
#                               'verticalalignment': 'bottom'
#                               })

#             ax.text((s+e)/2, h_s+h+.2, text,
#                     **kvargs.setdefault('kv_line_text',_fontdict))

#     ppara_plot_function_custom['heatmap_gene'] = _func

# with Block("ppara_plot_function_custom['heatmap_gene_get_marker_and_df_plot']"):
#     def _func(adata, group_by, marker, layer='scaled'):
#         marker = pd.concat([pd.DataFrame({
#             'cell_type': k,
#             'gene': list(vs)
#         }) for k, vs in marker.items()])
#         marker.index = np.arange(marker.shape[0])
#         marker['line_start'] = marker.groupby('cell_type').cumcount() == 0
#         marker['line_end'] = [False] + marker['line_start'].to_list()[1:]
#         df_plot = sc.get.obs_df(adata,[group_by] + list(marker['gene'].unique()),
#             layer=layer)
#         df_plot = df_plot.loc[:,df_plot.columns.unique()]
#         df_plot = group_agg(
#             df_plot, [group_by], {
#                 g: 'mean' for g in marker['gene']}, reindex=False)
#         df_plot = df_plot.loc[:,marker['gene']]
#         df_plot.index.name = ''
#         return marker, df_plot

#     ppara_plot_function_custom['heatmap_gene_get_marker_and_df_plot'] = _func

# with Block("ppara_plot_function_custom['heatmap_gene_process_multi_marker_and_df_plot']"):
#     def _func(list_data):
#         df_plot = [_[1].assign(**{'gap_{}'.format(_i): np.nan})
#                    for _i, (_) in enumerate(list_data)]
#         df_plot = pd.concat(df_plot, axis=1).copy()

#         df_temp = df_plot.columns.to_frame(name='gene')
#         df_temp['group'] = df_temp['gene'].str.extract('gap_(\d+)',expand=False)
#         df_temp['group'] = df_temp['group'].bfill()
#         df_temp.dtypes

#         df_marker = [_[0].assign(**{'group':str(_i)}) for _i, (_) in enumerate(list_data)]
#         df_marker = pd.concat(df_marker, axis=0)

#         # 防止两部分中存在同名的基因,添加group作为merge的条件
#         df_marker = pd.merge(df_temp,df_marker,
#             on='gene,group'.split(','),how='left')
#         df_marker['line_start'] = df_marker['line_start'].mask(
#             df_marker['line_start'].isna(), False)
#         df_marker['line_end'] = df_marker['line_end'].mask(
#             df_marker['gene'].str.match('gap_\\d+'), True
#         )

#         return df_marker, df_plot

#     ppara_plot_function_custom['heatmap_gene_process_multi_marker_and_df_plot'] = _func


# In[ ]:


with Block("ppara_plot_function_custom['heatmap_gene_get_marker_and_df_plot']"):
    def _func(adata, group_by, marker, layer='scaled'):
        marker = pd.concat([pd.DataFrame({
            'cell_type': k,
            'gene': list(vs)
        }) for k, vs in marker.items()])
        marker.index = np.arange(marker.shape[0])
        marker['line_start'] = marker.groupby('cell_type').cumcount() == 0
        marker['line_end'] = [False] + marker['line_start'].to_list()[1:]
        df_plot = sc.get.obs_df(
            adata,
            [group_by] +
            list(
                marker['gene'].unique()),
            layer=layer)
        df_plot = group_agg(
            df_plot, [group_by], {
                g: 'mean' for g in marker['gene']}, reindex=False)
        df_plot = df_plot.loc[:, marker['gene']]
        df_plot.columns = pd.MultiIndex.from_frame(
            marker.loc[:, 'cell_type,gene'.split(',')])
        return marker, df_plot

    ppara_plot_function_custom['heatmap_gene_get_marker_and_df_plot'] = _func

with Block("ppara_plot_function_custom['heatmap_gene_process_multi_marker_and_df_plot']"):
    def _func(list_data):
        def add_col_gap(df, key, values):
            df[key] = values
            return df

        df_plot = pd.concat([add_col_gap(_[1], ('gap', 'gap_{}'.format(
            _i)), np.nan) for _i, (_) in enumerate(list_data)], axis=1).copy()

        df_temp = df_plot.columns.to_frame()
        df_temp.index = np.arange(df_temp.shape[0])
        df_temp['group'] = df_temp['gene'].str.extract(
            'gap_(\\d+)', expand=False)
        df_temp['group'] = df_temp['group'].bfill()

        df_marker = df_marker = [_[0].assign(**{'group': str(_i)})
                                 for _i, (_) in enumerate(list_data)]
        df_marker = pd.concat(df_marker, axis=0)
        df_marker = pd.merge(
            df_temp,
            df_marker,
            on='group,cell_type,gene'.split(','),
            how='left')
        df_marker['line_start'] = df_marker['line_start'].mask(
            df_marker['line_start'].isna(), False)
        df_marker['line_end'] = df_marker['line_end'].mask(
            df_marker['gene'].str.match('gap_\\d+'), True)
        df_marker['cell_type_mask'] = df_marker['cell_type'].mask(
            df_marker['cell_type'].str.match('^gap$') &
            df_marker['gene'].str.match('^gap_\\d+$'), '')
        df_marker['gene_mask'] = df_marker['gene'].mask(
            df_marker['cell_type'].str.match('^gap$') &
            df_marker['gene'].str.match('^gap_\\d+$'), '')
        return df_marker, df_plot

    ppara_plot_function_custom['heatmap_gene_process_multi_marker_and_df_plot'] = _func

with Block("ppara_plot_function_custom['heatmap_gene']"):
    def _func(
            df_plot,
            df_marker,
            ax,
            zscore=True,
            cmap=cm,
            cbar=False,
            line_h=.25,
            kv_line={
                'linestyle': '-',
                'linewidth': 1,
                'color': 'grey',
                'dash_joinstyle': 'miter',
                'dash_capstyle': 'butt'},
            kv_heatmap={
                'vmax': 2.5,
                'vmin': -2.5},
            **kvargs):
        """
           kvargs:
               del_yticks: default False
               kv_line_text: default fontdict_default

"""
        df_marker.index = np.arange(df_marker.shape[0])
        df_temp = df_plot.columns.to_frame().copy()
        df_temp.index = np.arange(df_temp.shape[0])
        assert (df_marker.loc[:, 'cell_type,gene'.split(',')] == df_temp).all(
        ).all(), """[Error] df_plot.columns is not equal to
 df_marker.loc[:,'cell_type,gene'.split(',')]"""

        # if df_marker['gene'][-1].startswith('gap_'):
        if not df_marker['gene_mask'].to_numpy()[-1]:
            df_plot = df_plot.iloc[:, :-1].copy()
        if zscore:
            df_plot = stats.zscore(df_plot, axis=0, nan_policy='omit')
        sns.heatmap(df_plot, square=True, cbar=cbar,
                    cmap=cmap, ax=ax, **kv_heatmap)
        ax.set_xlabel(''), ax.set_ylabel('')
        if kvargs.setdefault('del_yticks', False):
            ax.set_yticks([], [], rotation=0,
                          **fontdict_default)
        else:
            _texts = df_plot.index
            ax.set_yticks(
                np.arange(
                    df_plot.shape[0])+.5,
                df_plot.index,
                rotation=0,
                **fontdict_default)
        _texts = df_marker[df_marker['gene_mask'].str.len()
                           > 0]['gene_mask']
        ax.set_xticks(_texts.index + .5, _texts, rotation=90,
                      **fontdict_default)
        ax.set_ylim(min(ax.get_ylim()), max(ax.get_ylim())+1)

        # plot line
        for text, s, e in zip(
                df_marker.query('line_start')['cell_type'],
                df_marker.query('line_start').index,
                df_marker.query('line_end').index.to_list() + [df_marker.index.size]):
            h_s, h = df_plot.shape[0] + .25, line_h
            s = s+.25
            e = e-.25
            # text = 'mast cell'
            ax.step([s, s, e, e], [h_s, h_s+h, h_s+h, h_s], **kv_line)
            _fontdict = fontdict_default.copy()
            _fontdict.update({'horizontalalignment': 'center',
                              'multialignment': 'center',
                              'verticalalignment': 'bottom'
                              })

            ax.text((s+e)/2, h_s+h+.2, text,
                    **kvargs.setdefault('kv_line_text', _fontdict))

    ppara_plot_function_custom['heatmap_gene'] = _func


# # 多重检验 multiple_test

# In[10]:


def multiple_test(data, key_groupby, key_value, test_pairs,
                  test_func, test_func_kwarg={}, fd_method='bh'):
    """多重检验
    """
    def _get_values(data, key_groupby, value, key_value):
        return data.query(
            "{} == '{}'".format(
                key_groupby,
                value))[key_value].to_numpy()

    def _group_false_discovery_control(
            df, key_groupby='value_x', fd_method='bh'):

        df_list = []
        for g in df[key_groupby].unique():
            pass

            _temp = res.query("{} == '{}'".format(key_groupby, g)).copy()
            _temp['padj'] = stats.false_discovery_control(
                _temp['pvalue'], method=fd_method)
            df_list.append(_temp)
        return pd.concat(df_list)

    # 处理 test_pairs
    if isinstance(test_pairs[0], str):
        test_pairs[0] = [test_pairs[0]]
    if len(test_pairs[0]) == 1:
        test_pairs[0] = [test_pairs[0][0]
                         for _ in range(len(test_pairs[1]))]
    test_pairs[0] = list(test_pairs[0])
    test_pairs[1] = list(test_pairs[1])
    _temp = np.setdiff1d(
        np.unique(
            np.ravel(test_pairs)),
        data[key_groupby].unique())
    assert _temp.size == 0, '[Error] not all element of test_pairs are in data[key_gourpby]\n\t{}'\
        .format(','.join(_temp))
    # 多重假设检验并校正pvalue
    res = pd.MultiIndex.from_arrays(test_pairs).to_frame(
        name='value_x,value_y'.split(','))
    res = res.drop_duplicates(
        'value_x,value_y'.split(',')).query("value_x != value_y")

    res['mean_x'] = res.apply(lambda row: np.mean(_get_values(
        data, key_groupby, row['value_x'], key_value)), axis=1)
    res['mean_y'] = res.apply(lambda row: np.mean(_get_values(
        data, key_groupby, row['value_y'], key_value)), axis=1)
    res['mean_diff'] = res.eval('mean_y - mean_x')
    res['percent_mean_diff'] = res.eval('mean_diff/mean_x * 100')

    res['test_res'] = res.apply(
        lambda row: test_func(
            x=_get_values(
                data,
                key_groupby,
                row['value_x'],
                key_value),
            y=_get_values(
                data,
                key_groupby,
                row['value_y'],
                key_value),
            **test_func_kwarg),
        axis=1)
    res['statistic'] = res['test_res'].apply(lambda x: x.statistic)
    res['pvalue'] = res['test_res'].apply(lambda x: x.pvalue)
    # 以value_x分组对pvalue进行校正
    res = _group_false_discovery_control(res, 'value_x', fd_method)
    return res.drop(columns=['test_res'])


# # other func

# In[11]:


def show_dict_key(data, tag='', sort_key=True):
    print("\n>{}['']".format(tag).ljust(75, '-'))
    ks = list(data.keys())
    if sort_key:
        ks = np.sort(ks)
    print(*['\t{}'.format(k) for k in ks], sep='\n')


def str_next_line(s, s_find, index=30):
    """将字符串s每index个字符分为一组,
除了第0组外，其余组将第一个s_find替换为\n
    """
    _left = np.arange(0, len(s), index)
    _right = np.concatenate([_left[1:], [len(s)]])
    return ''.join([s[_l:_r] if _i == 0 else s[_l:_r].replace(
        s_find, '\n', 1) for _i, (_l, _r) in enumerate(zip(_left, _right))])


def mpl_subplots_get_fig_axs(nrows=1, ncols=1,
                             ratio_nrows=2, ratio_ncols=2):
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(
        ratio_nrows*ncols, ratio_ncols*nrows)
    )
    axs = np.ravel(axs)
    if axs.size == 1:
        axs = axs[0]

    return fig, axs


def yield_coordinate_from_arr(arr, new_shape, order='c', show=False):
    # arr = np.arange(12)
    if show:
        display(np.reshape(arr, newshape=new_shape, order=order))
    x, y = np.where(
        np.reshape(
            arr, newshape=new_shape, order=order) > -1)  # 坐标
    for i in range(arr.size):
        yield (x[i], y[i], arr[i])

