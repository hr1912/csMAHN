#!/usr/bin/env python
# coding: utf-8

# ```shell
# conda activate
# cd ~/link/res_publish
# 
# jupyter nbconvert func.ipynb --to python
# 
# jupyter nbconvert func_r_map_seruat.ipynb --to python
# mv func_r_map_seruat.py func_r_map_seruat.r
# 
# jupyter nbconvert README.ipynb --to markdown
# 
# echo 'finish'
# 
# ```

# In[1]:


from pathlib import Path
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
from IPython.display import display

from zipfile import ZipFile
import tarfile
from collections import namedtuple
import time
import re
import warnings
from json import loads, dumps
from io import StringIO


import came
import csMAHN
from csMAHN.utils import utility  # 避免pp 时的循环导入
from csMAHN.utils import preprocess as pp
from csMAHN.utils.train import Trainer


# ## SAMap package and a new function

# In[2]:


import samap
from samalg import SAM
from samap.mapping import SAMAP
from samap import analysis as sana
import scipy as sp

map_sp_SAMap = {'human': 'hu',
                'mouse': 'mm',
                'zebrafish': 'zf',
                'chicken': 'ch',
                'macaque': 'ma'}


def get_alignment_score_for_each_cell(sm, keys, n_top=0):
    "n_top 无效 但保留"
    def customize_compute_csim(samap, key, X=None, prepend=True, n_top=0):
        splabels = sana.q(samap.adata.obs['species'])
        skeys = splabels[np.sort(
            np.unique(splabels, return_index=True)[1])]

        cl = []
        clu = []
        for sid in skeys:
            if prepend:
                cl.append(
                    sid +
                    '_' +
                    sana.q(
                        samap.adata.obs[key])[
                        samap.adata.obs['species'] == sid].astype('str').astype('object'))
            else:
                cl.append(sana.q(samap.adata.obs[key])[
                          samap.adata.obs['species'] == sid])
            clu.append(np.unique(cl[-1]))

        clu = np.concatenate(clu)
        cl = np.concatenate(cl)

        CSIM = np.zeros((clu.size, clu.size))
        if X is None:
            X = samap.adata.obsp["connectivities"].copy()

        xi, yi = X.nonzero()
        spxi = splabels[xi]
        spyi = splabels[yi]

        filt = spxi != spyi
        di = X.data[filt]
        xi = xi[filt]
        yi = yi[filt]

        px, py = xi, cl[yi]
        p = px.astype('str').astype('object')+';'+py.astype('object')

        A = pd.DataFrame(data=np.vstack((p, di)).T, columns=["x", "y"])
        valdict = sana.df_to_dict(A, key_key="x", val_key="y")
        cell_scores = [valdict[k].sum() for k in valdict.keys()]
        ixer = pd.Series(data=np.arange(clu.size), index=clu)
        if len(valdict.keys()) > 0:
            xc, yc = sana.substr(list(valdict.keys()), ';')
            xc = xc.astype('int')
            yc = ixer[yc].values
            cell_cluster_scores = sp.sparse.coo_matrix(
                (cell_scores, (xc, yc)), shape=(X.shape[0], clu.size)).A
            return pd.DataFrame(
                cell_cluster_scores,
                index=samap.adata.obs.index,
                columns=clu)
        else:
            raise Exception('[Error]')
            # return np.zeros((clu.size, clu.size)), clu

    if len(list(keys.keys())) < len(list(sm.sams.keys())):
        samap = SAM(counts=sm.samap.adata[np.in1d(
            sm.samap.adata.obs['species'], list(keys.keys()))])
    else:
        samap = sm.samap

    clusters = []
    ix = np.unique(samap.adata.obs['species'], return_index=True)[1]
    skeys = sana.q(samap.adata.obs['species'])[np.sort(ix)]

    for sid in skeys:
        clusters.append(sana.q([sid+'_'+str(x)
                        for x in sm.sams[sid].adata.obs[keys[sid]]]))

    cl = np.concatenate(clusters)
    l = "{}_mapping_scores".format(';'.join([keys[sid] for sid in skeys]))
    samap.adata.obs[l] = pd.Categorical(cl)

    # CSIMth, clu = _compute_csim(samap, l, n_top = n_top, prepend = False)
    cell_cluster_scores = customize_compute_csim(
        samap, l, n_top=n_top, prepend=False)
    return cell_cluster_scores
# cell_cluster_scores = get_alignment_score_for_each_cell(sm, keys, n_top = 0)
# cell_cluster_scores


# In[3]:


p_link = Path("/public/workspace/licanchengup/link")
p_publish = p_link.joinpath("res_publish")
p_run = p_publish.joinpath("run")
p_plot = p_publish.joinpath("plot")
p_res = p_publish.joinpath("res")
p_cache = p_run.joinpath("cache")
p_df_varmap = p_publish.joinpath('homo/df_varmap.csv')

map_sp = {k: v for k, v in zip(
    'h,m,z,ma,c'.split(','),
    'human,mouse,zebrafish,macaque,chicken'.split(',')
)}
map_sp_reverse = {v: k for k, v in map_sp.items()}



# # F1 score

# In[4]:


def calculate_confusion_matrix(y_true, y_pred):
    """
    计算混淆矩阵
    """
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    classes = np.sort(np.union1d(np.unique(y_pred), np.unique(y_true)))
    matrix = pd.DataFrame(0, index=classes, columns=classes)
    matrix.update(pd.DataFrame({
        'y_true': y_true,
        'y_pred': y_pred}
    ).groupby(['y_true', 'y_pred'])['y_pred'].count().to_frame('count').reset_index().pivot(
        index='y_true',
        columns='y_pred',
        values='count'
    ).fillna(0)
    )
    matrix.index.name = 'true'
    matrix.columns.name = 'pred'
    return matrix


def calculate_more_with_confusion_matrix(data):
    df = pd.DataFrame({'TP': np.diag(data),
                       'FP': data.sum(axis=0) - np.diag(data),
                       'FN': data.sum(axis=1) - np.diag(data)})
    df['Precision'] = df.eval('TP /(TP + FP)').fillna(0)
    df['Recall'] = df.eval('TP /(TP + FN)').fillna(0)
    df['F1 Score'] = df.eval(
        '2*(Precision * Recall) / (Precision + Recall)').fillna(0)
    return df


def calculate_accuracy_with_confusion_matrix(data):
    data = data.loc[:, data.columns.isin(data.index)]
    return np.diag(data).sum() / np.sum(data.values)


def calculate_F1Score_with_confusion_matrix(data, average='weighted'):
    """
    average:
        macro (default): 算数均值
        weighted : 加权均值，以实际为真值的数量(TP+FN) 为权重
        micro: ？？？？Sum statistics over all labels
            分别对各类的TP , FP, FN 求和,在计算TP / (1/2(FP + FN))
            咦,这个不需要分别计算各类的F1-score
    """

    res = None
    if average == 'macro':
        res = data['F1 Score'].mean()
    elif average == 'weighted':
        res = np.average(data['F1 Score'], weights=data.eval('TP+FN'))
    elif average == 'micro':
        res = data['TP'].sum() / (data['TP'].sum() + 1/2 *
                                  (data['FP'].sum() + data['FN'].sum()))
    else:
        raise Exception('[Error] average= {}'.format(average))
    return res


def is_1v1(row):
    res = None
    if row['model'] == 'Seurat':
        res = True
    elif row['model'] in ['came', 'csMAHN']:
        if 'is_1v1=True' in row['resdir_tag']:
            res = True
        if 'is_1v1=False' in row['resdir_tag']:
            res = False
    return res


def get_significance_marker(p_value, markers={
    '**': 0.001,
    '**': 0.01,
    '*': 0.05
}, not_significance_marker='ns'):

    res = not_significance_marker
    markers = pd.Series(markers)
    markers = markers[markers > p_value]
    if markers.size > 0:
        res = markers[markers == markers.min()].index[0]
    return res


# # plot

# In[5]:


def savefig(fig,fig_name,p_plot=p_plot):
    fig.savefig(p_plot.joinpath(fig_name), transparent=True, dpi=200,bbox_inches ='tight')
    print('[out][plot] {} \n\tin {}'.format(fig_name,p_plot))
    
def calculate_res_stat(row, q, update=False):
    p_res_stats = row['dir'].joinpath('res_stats.josn')
    res_stats = {}
    if p_res_stats.exists() and (not update):
        res_stats = loads(p_res_stats.read_text())
    else:
        df_obs = get_res_obs(row).query(q)
        df_obs.index = df_obs.index.astype(str)
        cm = calculate_confusion_matrix(
            df_obs['true_label'], df_obs['pre_label'])
        acc = calculate_accuracy_with_confusion_matrix(cm)
        f1 = calculate_F1Score_with_confusion_matrix(
            calculate_more_with_confusion_matrix(cm))
        res_stats.update({
            q: {'confusion_matrix': cm.to_json(orient='columns'),
                'confusion_matrix_more': calculate_more_with_confusion_matrix(cm).to_json(orient='columns'),
                'Accuracy': calculate_accuracy_with_confusion_matrix(cm),
                'F1-score': calculate_F1Score_with_confusion_matrix(
                calculate_more_with_confusion_matrix(cm))
                }
        })
        p_res_stats.write_text(dumps(res_stats))

    for k, v in res_stats.items():
        res_stats[k]['confusion_matrix'] = pd.read_json(
            StringIO(v['confusion_matrix']), orient='columns')
        res_stats[k]['confusion_matrix_more'] = pd.read_json(
            StringIO(v['confusion_matrix_more']), orient='columns')
    return res_stats[q]


def get_res_stat(row, q, key, update=False):
    res_stat = calculate_res_stat(row, q, update)
    return res_stat.setdefault(key, None)


def time_tag_detect(p):
    return True if re.match('.+;\\d{6}-\\d{4}$', Path(p).name) else False


def time_tag_toggle(p):
    p = Path(p)
    assert p.exists(), '[not exists]\n{}'.format(p)
    p_res = p
    if time_tag_detect(p):
        # raise Exception('[time tag has existsed]\n{}'.format(p.name))
        p_res = p.with_name(re.sub(';\\d{6}-\\d{4}$', '', p_res.name))
    else:
        p_res = p.with_name(
            '{};{}'.format(
                p.name,
                time.strftime(
                    '%y%m%d-%H%M',
                    time.localtime())))
    assert not p_res.exists(), '[target has existed]\n{}'.format(p_res)
    p.rename(p_res)
    return p_res


# # other functions

# In[ ]:


def group_agg(
        obs,
        groupby_list,
        agg_dict,
        dropna=True,
        reindex=True,
        rename_dict=None):
    res = obs.groupby(groupby_list, dropna=dropna).agg(agg_dict)
    if reindex:
        res.columns = ["_".join(i) for i in res.columns]
        res = res.index.to_frame().join(res)
        res.index = np.arange(res.shape[0])
    if isinstance(rename_dict, dict):
        res = res.rename(columns=lambda k: rename_dict.setdefault(k, k))
    return res


def rm_rf(p):
    if not p.exists():
        return

    if p.is_file():
        p.unlink()

    if p.is_dir():
        for i in p.iterdir():
            if i.is_file():
                i.unlink()
            if i.is_dir():
                rm_rf(i)  # 递归
        p.rmdir()


def h5ad_to_mtx(adata, p_dir, prefixes="", as_int=True):
    """
    将adata对象保存为mtx
    p_dir ： 输出路径
    as_int : 是否将矩阵转换int类型
        default True
    """
    assert adata.obs.index.is_unique, '[Error] obs index is not unique'
    assert adata.var.index.is_unique, '[Error] var index is not unique'

    p_dir = Path(p_dir)
    p_dir.mkdir(parents=True, exist_ok=True)

    # [out] genes.tsv
    adata.var["gene_names"] = adata.var_names.to_numpy()
    if "gene_ids" not in adata.var_keys():
        adata.var["gene_ids"] = adata.var["gene_names"]
    df_genes = adata.var.loc[:, ["gene_ids", "gene_names"]]
    df_genes.to_csv(
        p_dir.joinpath("{}genes.tsv".format(prefixes)),
        header=False,
        index=False,
        sep="\t",
    )

    # [out] barcodes.tsv obs.csv
    adata.obs.loc[:, []].to_csv(
        p_dir.joinpath("{}barcodes.tsv".format(prefixes)),
        header=False,
        index=True,
        sep="\t",
    )

    if len(adata.obs_keys()) > 0:
        adata.obs.to_csv(
            p_dir.joinpath("{}obs.csv".format(prefixes)), index=True
        )

    # [out] matrix.mtx
    adata.X = csr_matrix(adata.X)
    if as_int:
        adata.X = adata.X.astype(int)
    nonzero_index = [i[:10] for i in adata.X.nonzero()]
    print(
        "frist 10 data.X nonzero elements:\n",
        adata.X[nonzero_index[0], nonzero_index[1]],
    )
    mmwrite(
        p_dir.joinpath("{}matrix.mtx".format(prefixes)), adata.X.getH()
    )
    print("[out] {}".format(p_dir))

def get_path_varmap(
        sp_ref,
        sp_que,
        p_df_varmap=p_df_varmap,
        model='csMAHN'):
    """通过sp_ref 和 sp_que获取path_varmap
    ./homo/df_varmap.csv 存储了
    path_varmap路径及信息
"""
    if model in 'csMAHN,came'.split(','):
        if isinstance(p_df_varmap, str) or isinstance(p_df_varmap, Path):
            p_df_varmap = pd.read_csv(p_df_varmap)
        index_ = p_df_varmap.query("sp_ref == '{}' & sp_que == '{}'".format(
            sp_ref, sp_que)).index
        if index_.size == 1:
            res = Path(p_df_varmap.loc[index_[0], 'path'])
            if res.exists():
                return res
            else:
                raise Exception("[not exists] {}".format(res))
        else:
            raise Exception(
                "[get {} path]can not get speicifed and unique path\nsp_ref\tsp_que\n{}\t{}".format(
                    index_.size, sp_ref, sp_que))
    elif model == 'SAMap':
        return Path(p_df_varmap).parent.joinpath('SAMap/maps_gene_name')
    else:
        raise Exception(
            "[Error] can not find path_varmap with model '{}'".format(model))

def get_test_result_df(
        p,
        extract="^(?P<tissue>.+)_(?P<sp_ref>.+)-corss-(?P<sp_que>.+);(?P<model>came|csMAHN|Seurat|SAMap);(?P<name_ref>[\\w-]+)-map-(?P<name_que>[[\\w-]+);?(?P<resdir_tag>.+)?$"):
    p = Path(p)
    df = pd.DataFrame({"dir": [i for i in p.iterdir() if i.is_dir()]})
    df = df[df["dir"].apply(lambda x: x.joinpath("finish").exists())]
    df["name"] = df["dir"].apply(lambda x: x.name)
    print('\n[extract]\n{}'.format(extract))
    df = df.join(
        df["name"].str.extract(
            extract
        )
    )
    return df

def get_res_obs(row):
    """
    row 为来自get_test_result_df的行
    dataset,cell_type,true_label,pre_label,max_prob,is_right,UMAP1,UMAP2
    """
    def get_res_obs_csMAHN(row):
        df = pd.read_csv(
            row["dir"].joinpath("res_2", "pre_out_2.csv"), index_col=0
        ).loc[:, "true_label,pre_label,pre_label_NUM,max_prob".split(",")]

        p_am = row["dir"].joinpath("adata_meta")
        if not p_am.exists():
            print("[not exists] {}/{} ".format(row["dir"].name, p_am.name))
            # return df
        df_obs = pd.read_csv(
            row["dir"].joinpath("adata_meta", "obs.csv"), index_col=0
        ).loc[:, ["dataset", "cell_type"]]

        df_umap = pd.read_csv(
            p_am.joinpath("obsm.csv"), index_col=None
        ).rename(columns={"X_umap1": "UMAP1", "X_umap2": "UMAP2"})
        df_umap.index = df_obs.index

        if not df_obs.index.is_unique:
            if df_obs.apply(
                    lambda row: "{};{}".format(
                        row.name,
                        row['dataset']),
                    axis=1).is_unique:
                df_obs.index = df_obs.apply(
                    lambda row: "{};{}".format(
                        row.name, row['dataset']), axis=1).to_numpy()
                df.index = df_obs.index
                df_umap.index = df_obs.index
                print('[index] add dataset')
            else:
                raise Exception('[index not unique]')
        df = df.join(df_obs).join(df_umap)

        df["is_right"] = df.eval("true_label == pre_label")
        df = df.loc[
            :,
            np.intersect1d(
                "dataset,cell_type,true_label,pre_label,max_prob,is_right,UMAP1,UMAP2".split(
                    ","
                ),
                df.columns,
            ),
        ]
        return df

    def get_res_obs_came(row):
        df = pd.read_csv(
            row["dir"].joinpath("obs.csv"),
            index_col=0,
            usecols="original_name,dataset,REF,celltype,predicted,max_probs,is_right".split(
                ","
            ),
        )

        df.index = df.index.to_numpy()
        df = df.rename(
            columns={
                "REF": "true_label",
                "predicted": "pre_label",
                "max_probs": "max_prob",
            }
        ).rename(columns={"celltype": "cell_type"})

        df = df.loc[
            :,
            np.intersect1d(
                "dataset,cell_type,true_label,pre_label,max_prob,is_right,UMAP1,UMAP2".split(
                    ","
                ),
                df.columns,
            ),
        ]
        p_am = row["dir"].joinpath("adata_meta")
        if not p_am.exists():
            print("[not exists] {}/{} ".format(row["dir"].name, p_am.name))
            return df

        # get UMAP1,UMAP2
        temp_df = pd.read_csv(
            p_am.joinpath("obsm.csv"), index_col=None
        ).rename(columns={"X_umap1": "UMAP1", "X_umap2": "UMAP2"})
        assert (
            temp_df.shape[0] == df.shape[0]
        ), "[Error][length not equal] {} {}".format(
            temp_df.shape[0], df.shape[0]
        )
        temp_df.index = df.index
        df = df.join(temp_df)

        df["is_right"] = df.eval("true_label == pre_label")
        df = df.loc[
            :,
            np.intersect1d(
                "dataset,cell_type,true_label,pre_label,max_prob,is_right,UMAP1,UMAP2".split(
                    ","
                ),
                df.columns,
            ),
        ]
        return df

    def get_res_obs_Seurat(row):
        return pd.read_csv(
            row["dir"].joinpath("obs.csv"), index_col=0, low_memory=False
        ).loc[
            :,
            "UMAP1,UMAP2,cell_type,dataset,is_right,max_prob,pre_label,true_label".split(
                ","
            ),
        ]

    def get_res_obs_SAMap(row):
        return get_res_obs_Seurat(row)

    map_func = {
        k: v
        for k, v in zip('csMAHN,came,Seurat,SAMap'.split(','),
                        [get_res_obs_csMAHN, get_res_obs_came, get_res_obs_Seurat, get_res_obs_SAMap])
    }
    assert row["model"] in map_func.keys(
    ), "can note get res obs. model = {}".format(row["model"])

    res = map_func[row["model"]](row)
    res['dataset_type'] = res['dataset'].map({
        '{tissue}_{sp_ref}'.format(**row): 'ref',
        '{tissue}_{sp_que}'.format(**row): 'que'
    })
    return res

def show_umap(row):
    """
    row 为来自get_test_result_df的行

    show test result umap
    row
        umap_dataset      :  png path
        umap_umap         :  png path
        p_cell_type_table :  csv path
    """
    row['umap_dataset'] = Path(
        row["dir"]).joinpath(
        "figs",
        "umap_dataset.png")
    row['umap_umap'] = Path(row["dir"]).joinpath('figs', 'umap_umap.png')
    row['p_cell_type_table'] = Path(
        row["dir"]).joinpath('group_counts.csv')
    print(row['name'].ljust(75, '-'))
    fig, ax = plt.subplots(1, 2, figsize=(14, 7))
    ax[0].imshow(mpimg.imread(row["umap_dataset"]))
    ax[0].set_axis_off()
    ax[1].imshow(mpimg.imread(row["umap_umap"]))
    ax[1].set_axis_off()
    display(fig, pd.read_csv(row["p_cell_type_table"], index_col=0))
    fig.clear()


# # cross species Model

# In[ ]:


def load_adata(p):
    def load_h5ad_from_mtx(p):
        p = Path(p)
        assert p.joinpath("matrix.mtx").exists(
        ), '[not exists]matrix.mtx\nin{}'.format(p)
        adata = sc.read_10x_mtx(p)
        if p.joinpath("obs.csv").exists():
            adata.obs = pd.read_csv(p.joinpath("obs.csv"), index_col=0)
            adata.obs.to_numpy()
            # adata.obs = adata.obs.loc[:, []].join(
            #     pd.read_csv(p.joinpath("obs.csv"), index_col=0))
        else:
            print('[not exists]obs.csv\nin {}'.format(p))
        return adata
    p = Path(p)
    if p.match("*.h5ad"):
        return sc.read_h5ad(p)
    elif p.is_dir() and p.joinpath("matrix.mtx").exists():
        return load_h5ad_from_mtx(p)
    else:
        raise Exception("[can not load adata] {}".format(p))


def load_normalized_adata(p):
    p = Path(p)
    assert p.is_dir, '[Error] please get a path of dir'
    p_h5ad = p.joinpath('normalize.h5ad')
    adata = None
    if p_h5ad.exists():
        adata = sc.read_h5ad(p_h5ad)
    else:
        adata = load_adata(p)
        print('[normalize adata]\n{}\n'.format(p_h5ad))
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        adata.write_h5ad(p_h5ad)
    return adata


def get_1v1_matches(
        df_match,
        key_homology_type='homology_type',
        value_homology_type='ortholog_one2one'):
    """
    from came.pp.take_1v_matches
    """
    l, r = df_match.columns[:2]
    l_unique = df_match[l].value_counts(
    ).to_frame().query("count == 1").index
    r_unique = df_match[r].value_counts(
    ).to_frame().query("count == 1").index
    keep = pd.DataFrame({
        'l_is_unique': df_match[l].isin(l_unique),
        'r_is_unique': df_match[r].isin(r_unique)
    }).min(axis=1)
    df_match = df_match[keep]
    df_match = df_match.query(
        "{} == '{}'".format(
            key_homology_type,
            value_homology_type))
    return df_match


def get_homology_parameters(adata1, adata2, df_varmap):
    res = {}
    df_homo_paras = pd.DataFrame({
        "gn_ref": adata1.var_names
    }).merge(df_varmap, on='gn_ref', how='left')
    res['homology_one2one_find'] = df_homo_paras['gn_que'].notna().sum()

    df_homo_paras['gn_que'] = df_homo_paras.apply(
        lambda row: 'not_o2o_' + row['gn_ref'] if pd.isna(
            row['gn_que']) else row['gn_que'], axis=1)

    assert (df_homo_paras['gn_ref'] == adata1.var_names).all(
    ), "df_homo_paras['gn_ref'] not equal adata1.var_names"
    assert df_homo_paras['gn_ref'].is_unique & df_homo_paras['gn_que'].is_unique, "df_homo_paras gn_ref or gn_que is not unique"
    # came 和 csMAHN不做替换
    # adata1.var.index = df_homo_paras['gn_que'].to_numpy()
    res['homology_one2one_use'] = np.intersect1d(
        df_homo_paras['gn_que'],
        adata2.var.index).size
    res = {k: int(v) for k, v in res.items()}
    return res


def get_type_counts_info(adatas, key_class, dsnames):
    type_counts_list = []
    for i in range(len(adatas)):
        type_counts_list.append(pd.value_counts(adatas[i].obs[key_class]))
    counts_info = pd.concat(type_counts_list, axis=1, keys=dsnames)
    return counts_info


def aligned_type(adatas, key_calss):
    adata1 = adatas[0].copy()
    adata2 = adatas[1].copy()
    counts_info = get_type_counts_info(
        adatas, key_calss, dsnames=["reference", "query"]
    )
    print("----raw----")
    print(counts_info)
    counts_info = counts_info.dropna(how="any")
    print("----new----")
    print(counts_info)

    com_type = counts_info.index.tolist()
    adata1 = adata1[adata1.obs[key_calss].isin(com_type)]
    adata2 = adata2[adata2.obs[key_calss].isin(com_type)]
    return adata1, adata2


def unify_group_counts_index_name(resdir):
    """
    不知为何，group_counts.index.name参差不齐
    故将group_counts_unalign.index.name 赋给 group_counts.index.name
    """
    p_group_counts = resdir.joinpath("group_counts.csv")
    p_group_counts_unalign = resdir.joinpath("group_counts_unalign.csv")

    lines = p_group_counts.read_text().split("\n")
    lines[0] = ",".join(
        [p_group_counts_unalign.read_text().split("\n")[0].split(",")[0]]
        + lines[0].split(",")[1:]
    )
    p_group_counts.write_text("\n".join(lines))


# ## came

# In[ ]:


def precess_after_came(resdir, tissue_name, sp1, sp2, is_display=False):
    unify_group_counts_index_name(resdir)

    figdir = resdir / "figs"
    sc.settings.figdir = figdir

    display(
        pd.read_csv(resdir / "group_counts.csv", index_col=0)
    ) if is_display else None
    obs = pd.read_csv(resdir / "obs.csv", index_col=0)

    # umap
    # the last layer of hidden states
    h_dict = came.load_hidden_states(resdir / "hidden_list.h5")[-1]
    adt = came.pp.make_adata(
        h_dict["cell"], obs=obs, assparse=False, ignore_index=True
    )
    sc.pp.neighbors(adt, n_neighbors=15, metric="cosine", use_rep="X")
    sc.tl.umap(adt)

    sc.pl.umap(adt, color="dataset", save="_dataset.png")
    sc.pl.umap(adt, color="celltype", save="_umap.png")
    adt.write_csvs(resdir.joinpath("adata_meta"))
    # # umap.csv umap坐标存储
    # pd.DataFrame(
    #     adt.obsm["X_umap"],
    #     columns=["umap_1", "umap_2"],
    #     index=sc.get.obs_df(adt, "original_name"),
    # ).reset_index().to_csv(resdir.joinpath("umap.csv"), index=False)

    # test umap.csv
    # temp_obs = pd.read_csv(resdir.joinpath("obs.csv"))
    # temp_obs = temp_obs.merge(pd.read_csv(resdir.joinpath("umap.csv")),on="original_name")
    # temp_adata = sc.AnnData(obs=temp_obs)
    # sc.pl.scatter(temp_adata,"umap_1", "umap_2",color="dataset")
    # sc.pl.scatter(temp_adata,"umap_1", "umap_2",color="celltype")
    # del temp_adata,temp_obs

    # ratio
    display(
        obs["is_right"].sum() / obs["is_right"].size
    ) if is_display else None

    obs["name"] = obs["original_name"].str.extract(";(.+)", expand=False)
    obs["sp"] = obs["dataset"].str.extract(
        "%s_(\\w+)" % tissue_name, expand=False
    )
    # 物种
    df_sp = (
        group_agg(obs, ["sp"], {"is_right": ["sum", "count"]})
        .eval("ratio = is_right_sum/is_right_count")
        .sort_values(["sp"])
    )
    df_sp["type"] = "species"
    # 各个dataset
    df_dataset = (
        group_agg(obs, ["sp", "name"], {"is_right": ["sum", "count"]})
        .eval("ratio = is_right_sum/is_right_count")
        .sort_values(["sp", "name"])
    )
    df_dataset["type"] = "dataset"
    df_ratio = pd.concat([df_sp, df_dataset], axis=0)
    df_ratio["tissue"] = tissue_name
    df_ratio = df_ratio.loc[
        :,
        "tissue,type,sp,name,is_right_sum,is_right_count,ratio".split(","),
    ]
    display(df_ratio) if is_display else None
    df_ratio.to_csv(resdir / "ratio.csv", index=False)
    del df_sp, df_dataset, df_ratio

    # heatmap
    # all,reference,query
    for q, sp in zip(
        (
            "{col} == {col}".format(col=obs.columns[0]),
            # all, both sp1 and sp2
            "sp == '{sp}'".format(sp=sp1),
            "sp == '{sp}'".format(sp=sp2),
        ),
        ("all", sp1, sp2),
    ):
        res = (
            group_agg(
                obs.query(q),
                ["celltype", "predicted"],
                {"predicted": ["count"]},
            )
            .pivot(
                index="celltype",
                columns="predicted",
                values="predicted_count",
            )
            .fillna(0)
        )
        display(res) if is_display else None
        res.to_csv(resdir / f"predicted_count_{sp}.csv", index=True)
        ax = sns.heatmap(
            data=stats.zscore(res, axis=1), cmap=plt.get_cmap("Greens")
        )
        ax.set_title(f"{tissue_name}-{sp}")
        ax.figure.savefig(
            figdir / f"heatmap_ratio_{ax.get_title()}.pdf",
            bbox_inches="tight",
            dpi=120,
        )
        ax.figure.clear()


def run_came(
    path_adata1,
    path_adata2,
    key_class1,
    key_class2,
    sp1,
    sp2,
    tissue_name,
    path_varmap,
    aligned=False,
    resdir_tag=".",
    resdir=Path(
        "/public/workspace/licanchengup/download/test/test_result"
    ),
    limite_func=lambda adata1, adata2: (adata1, adata2), **kvargs
):
    """
    version:0.0.5
    kvargs:
        n_epochs: int
            default,500,但是为了与csMAHN统一
            n_epochs = sum(kvargs.setdefault("n_epochs",[100,200,200]))

        is_1v1: bool
            default,False

    """

    # Parameter settings
    # n_epochs = 500
    n_epochs = sum(kvargs.setdefault("n_epochs", [100, 200, 200]))
    batch_size = None
    n_pass = 100
    use_scnets = True
    ntop_deg = 50
    ntop_deg_nodes = 50
    node_source = "deg,hvg"
    # keep_non1v1_feats = True
    keep_non1v1_feats = not kvargs.setdefault("is_1v1", False)

    # setting directory for results
    if len(resdir_tag) > 0:
        resdir_tag = f"{tissue_name}_{sp1}-corss-{sp2};{resdir_tag}"
    else:
        resdir_tag = f"{tissue_name}_{sp1}-corss-{sp2}"

    resdir = resdir / resdir_tag

    # 终止 判断
    p_finish = resdir.joinpath("finish")
    if p_finish.exists():
        # precess_after_came(resdir,tissue_name,sp1, sp2)
        print(
            "[has finish]{} {}".format(
                time.strftime('%y%m%d-%H%M', time.localtime()),
                resdir.name)
        )
        return
    print(
        "[start]{} {}".format(
            time.strftime('%y%m%d-%H%M', time.localtime()),
            resdir.name

        ))
    # return

    figdir = resdir / "figs"
    sc.settings.figdir = figdir
    resdir.mkdir(parents=True, exist_ok=True)

    finish_content = ["[strat] %f" % time.time()]

    # # setting
    dsnames = (f"{tissue_name}_{sp1}", f"{tissue_name}_{sp2}")
    dsn1, dsn2 = dsnames
    homo_method = "biomart"

    # load data
    adata_raw1 = load_adata(path_adata1)
    adata_raw2 = load_adata(path_adata2)

    # limite 进一步对adata进行限制，默认不操作直接返回
    adata_raw1, adata_raw2 = limite_func(adata_raw1, adata_raw2)

    # group_counts_unalign.csv
    pd.concat(
        [
            adata_raw1.obs[key_class1].value_counts(),
            adata_raw2.obs[key_class2].value_counts(),
        ],
        axis=1,
        keys=dsnames,
    ).to_csv(resdir.joinpath("group_counts_unalign.csv"), index=True)
    # align
    if aligned:
        adata_raw1, adata_raw2 = aligned_type(
            [adata_raw1, adata_raw2], key_calss=key_class1
        )

    # 保存obs ,即真正测试的细胞的mata
    adata_raw1.obs.to_csv(resdir.joinpath("obs_ref.csv"), index=True)
    adata_raw2.obs.to_csv(resdir.joinpath("obs_que.csv"), index=True)

    # group_counts.csv
    temp = pd.concat(
        [
            adata_raw1.obs[key_class1].value_counts(),
            adata_raw2.obs[key_class2].value_counts(),
        ],
        axis=1,
        keys=dsnames,
    )
    # came会自行导出
    # temp.to_csv(resdir.joinpath("group_counts.csv"), index=True)
    # if temp.shape[0] < 2:
    #     # 错误标记
    #     print("[Error][group_counts no any item]")
    #     finish_content.append(
    #         "[Error][group_counts no any item] %f" % time.time()
    #     )
    #     p_finish.with_name("error").write_text("\n".join(finish_content))
    #     return

    adatas = [adata_raw1, adata_raw2]
    print("cell count --> %d" % sum([i.shape[0] for i in adatas]))

    df_varmap = pd.read_csv(path_varmap, usecols=range(3))
    df_varmap.columns = ["gn_ref", "gn_que", "homology_type"]
    if kvargs.setdefault("is_1v1", False):
        df_varmap = get_1v1_matches(df_varmap)
        homology_parameter = get_homology_parameters(
            adata_raw1, adata_raw2, df_varmap)
        print("""
[homology one2one]find {homology_one2one_find} genes
[homology one2one]use {homology_one2one_use} genes""".format(
            **homology_parameter))
        kvargs.update(homology_parameter)

    df_varmap_1v1 = came.pp.take_1v1_matches(df_varmap)

    kvargs.update({'path_adata1': str(path_adata1),
                   'path_adata2': str(path_adata2),
                   'key_class1': key_class1,
                   'key_class2': key_class2,
                   'sp1': sp1,
                   'sp2': sp2,
                   'tissue_name': tissue_name,
                   'path_varmap': str(path_varmap),
                   'aligned': aligned,
                   'resdir_tag': resdir_tag,
                   'resdir': str(resdir)})
    finish_content.append("[finish before run] %f" % time.time())
    warnings.filterwarnings("ignore")

    came_inputs, (adata1, adata2) = came.pipeline.preprocess_unaligned(
        adatas,
        key_class=key_class1,
        use_scnets=use_scnets,
        ntop_deg=ntop_deg,
        ntop_deg_nodes=ntop_deg_nodes,
        node_source=node_source,
    )

    outputs = came.pipeline.main_for_unaligned(
        **came_inputs,
        df_varmap=df_varmap,
        df_varmap_1v1=df_varmap_1v1,
        dataset_names=dsnames,
        key_class1=key_class1,
        key_class2=key_class2,
        do_normalize=True,
        keep_non1v1_feats=keep_non1v1_feats,
        n_epochs=n_epochs,
        resdir=resdir,
        n_pass=n_pass,
        batch_size=batch_size,
        plot_results=True,
    )

    finish_content.append("[finish run] %f" % time.time())

    dpair = outputs["dpair"]
    trainer = outputs["trainer"]
    h_dict = outputs["h_dict"]
    out_cell = outputs["out_cell"]
    predictor = outputs["predictor"]

    obs_ids1, obs_ids2 = dpair.obs_ids1, dpair.obs_ids2
    obs = dpair.obs
    classes = predictor.classes
    # 后处理
    precess_after_came(resdir, tissue_name, sp1, sp2)
    finish_content.append("[finish after run] %f" % time.time())

    # 完成标记
    resdir.joinpath("kvargs.json").write_text(dumps(kvargs))
    finish_content.append("[end] %f" % time.time())
    p_finish.write_text("\n".join(finish_content))


# ## csMAHN

# In[ ]:


def precess_after_csMAHN(
    resdir, tissue_name, sp1, sp2, is_display=False, **kvargs
):
    """
    kvargs:
        adt
    """
    unify_group_counts_index_name(resdir)

    assert resdir.joinpath("res_2").exists(), "[not exists] res_2"
    figdir = resdir / "figs"
    sc.settings.figdir = figdir

    display(
        pd.read_csv(resdir / "group_counts.csv", index_col=0)
    ) if is_display else None

    # 为obs添加dataset
    obs = pd.read_csv(
        resdir.joinpath("res_2", "pre_out_2.csv"), index_col=0
    )
    if not obs.index.is_unique:
        print("[obs index is not unique]")

    pre_obs = {
        f"{tissue_name}_{k}": pd.read_csv(
            resdir.joinpath(file_name), index_col=0
        )
        for k, file_name in zip([sp1, sp2], ["obs_ref.csv", "obs_que.csv"])
    }
    for k, value in pre_obs.items():
        value["dataset"] = k
        pre_obs[k] = value.loc[:, ["dataset"]]
    obs = pd.concat(pre_obs.values()).join(obs)
    obs["is_right"] = obs["true_label"] == obs["pre_label"]
    del pre_obs

    # umap
    adt = kvargs.setdefault("adt", None)
    assert adt is not None
    sc.pp.neighbors(adt, n_neighbors=15, metric="cosine", use_rep="X")
    sc.tl.umap(adt)
    sc.pl.umap(adt, color="dataset", save="_dataset.png")
    sc.pl.umap(adt, color="cell_type", save="_umap.png")

    adt.write_csvs(resdir.joinpath("adata_meta"))
    # # 存储adt.X 和 umap坐标
    # pd.DataFrame(adt.X, index=adt.obs.index).to_orc(resdir.joinpath("adt.X.orc"))
    # pd.DataFrame(
    #     adt.obsm["X_umap"], columns=["umap_1", "umap_2"], index=adt.obs.index
    # ).assign(**{k: adt.obs[k] for k in adt.obs.keys()}).reset_index().to_csv(
    #     resdir.joinpath("umap.csv"), index=False
    # )

    # # test umap.csv
    # temp_obs = pd.read_csv(resdir.joinpath("res_2", "pre_out_2.csv"), index_col=0)
    # temp_obs = temp_obs.join(pd.read_csv(resdir.joinpath("umap.csv"),index_col=0))
    # temp_adata = sc.AnnData(obs=temp_obs)
    # sc.pl.scatter(temp_adata,"umap_1", "umap_2",color="dataset")
    # sc.pl.scatter(temp_adata,"umap_1", "umap_2",color="cell_type")
    # del temp_adata,temp_obs

    # ratio
    obs["name"] = obs.index.to_series().str.extract(
        ";([^;]+)$", expand=False
    )
    obs["sp"] = obs["dataset"].str.extract(
        "%s_(\\w+)" % tissue_name, expand=False
    )
    display(
        obs["is_right"].sum() / obs["is_right"].size
    ) if is_display else None
    # 物种
    df_sp = group_agg(obs, ["sp"], {"is_right": ["count", "sum"]})
    df_sp["type"] = "species"
    # 各个dataset
    df_dataset = group_agg(
        obs, ["sp", "name"], {"is_right": ["count", "sum"]}
    )
    df_dataset["type"] = "dataset"

    df_ratio = pd.concat([df_sp, df_dataset]).eval(
        "ratio = is_right_sum/is_right_count"
    )
    df_ratio["tissue"] = tissue_name
    df_ratio = df_ratio.loc[
        :,
        "tissue,type,sp,name,is_right_sum,is_right_count,ratio".split(","),
    ]
    df_ratio.to_csv(resdir / "ratio.csv", index=False)
    del df_sp, df_dataset, df_ratio

    # heatmap
    # all,reference,query
    for q, sp in zip(
        (
            "{col} == {col}".format(col=obs.columns[0]),
            # all, both sp1 and sp2
            "sp == '{sp}'".format(sp=sp1),
            "sp == '{sp}'".format(sp=sp2),
        ),
        ("all", sp1, sp2),
    ):
        res = (
            group_agg(
                obs.query(q),
                ["true_label", "pre_label"],
                {"pre_label": ["count"]},
            )
            .pivot(
                index="true_label",
                columns="pre_label",
                values="pre_label_count",
            )
            .fillna(0)
        )
        display(res) if is_display else None
        res.to_csv(resdir / f"predicted_count_{sp}.csv", index=True)

        ax = sns.heatmap(
            data=stats.zscore(res, axis=1), cmap=plt.get_cmap("Greens")
        )
        ax.set_title(f"{tissue_name}-{sp}")
        ax.figure.savefig(
            figdir / f"heatmap_ratio_{ax.get_title()}.pdf",
            bbox_inches="tight",
            dpi=120,
        )
        ax.figure.clear()


def run_csMAHN(
    path_adata1,
    path_adata2,
    key_class1,
    key_class2,
    sp1,
    sp2,
    tissue_name,
    path_varmap,
    aligned=False,
    resdir_tag=".",
    resdir=Path('.'),
    limite_func=lambda adata1, adata2: (adata1, adata2),

    **kvargs
):
    """
    version = 0.0.9
    kvargs:
        n_epochs:
            default,[100, 200, 300]
            stages,即res_0,res_1，res_2 的 epochs
            累加制，res_0,res_1，res_2,实际epochs分别为100,300,600
            故最终epochs为stages之和
            stages = kvargs.setdefault("n_epochs",[100, 200, 300])


        is_1v1: bool
            default,False
    """
    homo_method = 'biomart'
    n_hvgs = 2000
    n_degs = 50
    seed = 123
    stages = kvargs.setdefault(
        'n_epochs', [
            100, 200, 200])  # [200, 200, 200]
    nfeats = 64  # enbedding size #128
    hidden = 64  # 128
    input_drop = 0.2
    att_drop = 0.2
    residual = True

    threshold = 0.9  # 0.8
    lr = 0.01  # lr = 0.01
    weight_decay = 0.001
    patience = 100
    enhance_gama = 10
    simi_gama = 0.1

    path_specie_1 = path_adata1
    path_specie_2 = path_adata2
    tissue = tissue_name
    species = [sp1, sp2]
    dsnames = (f"{tissue_name}_{sp1}", f"{tissue_name}_{sp2}")
    assert key_class1 == key_class2, "key_class is not equal"
    key_class = key_class1

    # make file to save
    resdir_tag = f"{tissue_name}_{sp1}-corss-{sp2};{resdir_tag}" if len(
        resdir_tag) > 0 else f"{tissue_name}_{sp1}-corss-{sp2}"
    curdir = os.path.join(resdir, resdir_tag)
    resdir = Path(curdir)
    model_dir = os.path.join(curdir, 'model_')
    figdir = os.path.join(curdir, 'figs')
    [Path(_).mkdir(exist_ok=True, parents=True)
     for _ in [curdir, model_dir, figdir]]
    for i in range(len(stages)):
        res_dir = os.path.join(curdir, f'res_{i}')
        Path(res_dir).mkdir(exist_ok=True, parents=True)
    checkpt_file = os.path.join(model_dir, "mutistages")
    print(checkpt_file)

    # is finish
    p_finish = Path(resdir).joinpath("finish")
    if p_finish.exists():
        print(
            "[has finish]{} {}".format(
                time.strftime('%y%m%d-%H%M', time.localtime()), resdir.name
            ))
        return
    print(
        "[start]{} {}".format(
            time.strftime('%y%m%d-%H%M', time.localtime()),
            resdir.name
        ))

    finish_content = ["[strat] %f" % time.time()]
    print('[path_varmap] {}'.format(path_varmap))
    adata_species_1 = load_adata(path_specie_1)
    adata_species_2 = load_adata(path_specie_2)

    # limite 进一步对adata进行限制，默认不操作直接返回
    adata_species_1, adata_species_2 = limite_func(
        adata_species_1, adata_species_2
    )
    # group_counts_unalign.csv
    pd.concat([adata_species_1.obs[key_class].value_counts(),
               adata_species_2.obs[key_class].value_counts(),],
              axis=1, keys=dsnames,).to_csv(
        resdir.joinpath("group_counts_unalign.csv"), index=True
    )
    # 仅保留公共细胞类群
    if aligned:
        adata_species_1, adata_species_2 = pp.aligned_type(
            [adata_species_1, adata_species_2], key_class
        )

    # group_counts.csv
    temp = pd.concat([adata_species_1.obs[key_class].value_counts(),
                      adata_species_2.obs[key_class].value_counts(),],
                     axis=1, keys=dsnames)
    print(temp)
    temp.to_csv(resdir.joinpath("group_counts.csv"), index=True)
    # if temp.shape[0] < 2:
    #     # 错误标记
    #     print("[Error][group_counts no any item]")
    #     finish_content.append(
    #         "[Error][group_counts no any item] %f" % time.time()
    #     )
    #     p_finish.with_name("error").write_text("\n".join(finish_content))
    #     return

    adata_species_1.obs.to_csv(resdir.joinpath("obs_ref.csv"), index=True)
    adata_species_2.obs.to_csv(resdir.joinpath("obs_que.csv"), index=True)

    # homo = pd.read_csv(path_varmap)
    homo = pd.read_csv(path_varmap, usecols=range(3))
    homo.columns = ["gn_ref", "gn_que", "homology_type"]
    if kvargs.setdefault("is_1v1", False):
        homo = get_1v1_matches(homo)
        homology_parameter = get_homology_parameters(
            adata_species_1, adata_species_2, homo)
        print("""
[homology one2one]find {homology_one2one_find} genes
[homology one2one]use {homology_one2one_use} genes""".format(
            **homology_parameter))
        kvargs.update(homology_parameter)

    kvargs.update({'path_adata1': str(path_adata1),
                   'path_adata2': str(path_adata2),
                   'key_class1': key_class1,
                   'key_class2': key_class2,
                   'sp1': sp1,
                   'sp2': sp2,
                   'tissue_name': tissue_name,
                   'path_varmap': str(path_varmap),
                   'aligned': aligned,
                   'resdir_tag': resdir_tag,
                   'resdir': str(resdir)})

    print(
        """Task: refernece:{} {} cells x {} gene -> query:{} {} cells x {} gene in {}""".format(
            dsnames[0],
            adata_species_1.shape[0],
            adata_species_1.shape[1],
            dsnames[1],
            adata_species_2.shape[0],
            adata_species_2.shape[1],
            tissue))

    start = time.time()
    finish_content.append("[finish before run] %f" % time.time())
    # knn时间较长
    print("\n[process_for_graph]\n".center(100, '-'))
    adatas, features_genes, nodes_genes, scnets, one2one, n2n = pp.process_for_graph(
        [adata_species_1, adata_species_2], homo, key_class, 'leiden', n_hvgs=n_hvgs, n_degs=n_degs)
    g, inter_net, one2one_gene_nodes_net, cell_label, n_classes, list_idx = pp.make_graph(adatas,
                                                                                          aligned,
                                                                                          key_class,
                                                                                          features_genes,
                                                                                          nodes_genes,
                                                                                          scnets,
                                                                                          one2one,
                                                                                          n2n,
                                                                                          has_mnn=True,
                                                                                          seed=seed)
    end = time.time()
    # 包括预处理时间
    print('Times preprocess for graph:{:.2f}'.format(end - start))
    print("\n[Trainer]\n".center(100, '-'))
    trainer = Trainer(adatas,
                      g,
                      inter_net,
                      list_idx,
                      cell_label,
                      n_classes,
                      threshold=threshold,
                      key_class=key_class)
    print("\n[train]\n".center(100, '-'))
    trainer.train(curdir=curdir,
                  checkpt_file=checkpt_file,
                  nfeats=nfeats,
                  hidden=hidden,
                  enhance_gama=enhance_gama,
                  simi_gama=simi_gama)

    finish_content.append("[finish run] %f" % time.time())
    adt = sc.AnnData(
        trainer.embedding_hidden.detach().numpy(),
        obs=pd.concat(
            [
                adatas[0]
                .obs.loc[:, [key_class]]
                .assign(dataset=dsnames[0]),
                adatas[1]
                .obs.loc[:, [key_class]]
                .assign(dataset=dsnames[1]),
            ]
        ).rename(columns={key_class: "cell_type"}),
    )
    # plot_umap(trainer.embedding_hidden, adatas, dsnames, figdir)
    precess_after_csMAHN(
        resdir, tissue_name, sp1, sp2, is_display=False, adt=adt
    )
    # 完成标记
    resdir.joinpath("kvargs.json").write_text(dumps(kvargs))
    finish_content.append("[end] %f" % time.time())
    p_finish.write_text("\n".join(finish_content))
    # return trainer, adatas


# ## SAMap

# In[ ]:


def run_SAMap(
    path_adata1,
    path_adata2,
    key_class1,
    key_class2,
    sp1,
    sp2,
    tissue_name,
    path_varmap,
    aligned=False,
    resdir_tag=".",
    resdir=Path('.'),
    limite_func=lambda adata1, adata2: (adata1, adata2),
    **kvargs
):

    path_specie_1 = path_adata1
    path_specie_2 = path_adata2
    tissue = tissue_name
    species = [sp1, sp2]
    dsnames = (f"{tissue_name}_{sp1}", f"{tissue_name}_{sp2}")
    assert key_class1 == key_class2, "key_class is not equal"
    key_class = key_class1

    # make file to save
    resdir_tag = f"{tissue_name}_{sp1}-corss-{sp2};{resdir_tag}" if len(
        resdir_tag) > 0 else f"{tissue_name}_{sp1}-corss-{sp2}"
    curdir = os.path.join(resdir, resdir_tag)
    resdir = Path(curdir)
    figdir = os.path.join(curdir, 'figs')
    [Path(_).mkdir(exist_ok=True, parents=True)
     for _ in [curdir, figdir]]

    # is finish
    p_finish = Path(resdir).joinpath("finish")
    if p_finish.exists():
        print(
            "[has finish]{} {}".format(
                time.strftime('%y%m%d-%H%M', time.localtime()),
                resdir.name))
        return
    print(
        "[start]{} {}".format(
            time.strftime('%y%m%d-%H%M', time.localtime()),
            resdir.name
        ))

    finish_content = ["[strat] %f" % time.time()]
    print('[path_varmap] {}'.format(path_varmap))
    adata_1 = load_adata(path_specie_1)
    adata_2 = load_adata(path_specie_2)
    assert pd.Series(
        np.concatenate(
            (adata_1.obs.index, adata_2.obs.index))).is_unique, '[Error] index is not unique'

    # limite 进一步对adata进行限制，默认不操作直接返回
    adata_1, adata_2 = limite_func(
        adata_1, adata_2
    )
    # group_counts_unalign.csv
    pd.concat([adata_1.obs[key_class].value_counts(),
               adata_2.obs[key_class].value_counts(),],
              axis=1, keys=dsnames,).to_csv(
        resdir.joinpath("group_counts_unalign.csv"), index=True
    )
    # 仅保留公共细胞类群
    if aligned:
        adata_1, adata_2 = pp.aligned_type(
            [adata_1, adata_2], key_class
        )

    # group_counts.csv
    temp = pd.concat([adata_1.obs[key_class].value_counts(),
                      adata_2.obs[key_class].value_counts(),],
                     axis=1, keys=dsnames)
    print(temp)
    temp.to_csv(resdir.joinpath("group_counts.csv"), index=True)
    # if temp.shape[0] < 2:
    #     # 错误标记
    #     print("[Error][group_counts no any item]")
    #     finish_content.append(
    #         "[Error][group_counts no any item] %f" % time.time()
    #     )
    #     p_finish.with_name("error").write_text("\n".join(finish_content))
    #     return

    adata_1.obs.to_csv(resdir.joinpath("obs_ref.csv"), index=True)
    adata_2.obs.to_csv(resdir.joinpath("obs_que.csv"), index=True)

    kvargs.update({'path_adata1': str(path_adata1),
                   'path_adata2': str(path_adata2),
                   'key_class1': key_class1,
                   'key_class2': key_class2,
                   'sp1': sp1,
                   'sp2': sp2,
                   'tissue_name': tissue_name,
                   'path_varmap': str(path_varmap),
                   'aligned': aligned,
                   'resdir_tag': resdir_tag,
                   'resdir': str(resdir)})

    start = time.time()
    finish_content.append("[finish before run] %f" % time.time())
    # SAMap -----------------------------------------------------

    sq_ref_SAMap = map_sp_SAMap[map_sp[sp1]]
    sq_que_SAMap = map_sp_SAMap[map_sp[sp2]]
    _SAMP = {
        sq_ref_SAMap: SAM(counts=adata_1),
        sq_que_SAMap: SAM(counts=adata_2)
    }
    [v.preprocess_data() for k, v in _SAMP.items()]
    [v.run() for k, v in _SAMP.items()]

    keys = {
        sq_ref_SAMap: key_class1,
        sq_que_SAMap: key_class2

    }
    sm = SAMAP(
        _SAMP,
        # f_maps + "{}{}/{}_to_{}.txt".format(id2, id1, id2, id1)
        # 。。。昂,不能传Path,得传str,还得加个 os.sep
        f_maps=str(path_varmap) if str(path_varmap).endswith(
            os.sep) else str(path_varmap) + os.sep,
        keys=keys

    )

    # run SAMap
    sm.run(pairwise=False)
    samap = sm.samap  # SAM object with three species stitched together

    finish_content.append("[finish run] %f" % time.time())

    alignment_score = samap.adata.obs.loc[:, ['species']].join(
        get_alignment_score_for_each_cell(sm, keys))
    alignment_score.to_csv(resdir.joinpath(
        'alignment_score_for_each_cell.csv'), index=True)
    alignment_score.head(2)

    alignment_score_query = alignment_score\
        .query("species == '{}'".format(sq_que_SAMap))\
        .filter(regex="^{}".format(sq_ref_SAMap))
    # 取每个细胞 在各类型上alignment score 最大值对应的类型
    alignment_score_query['pre_label'] = [
        alignment_score_query.columns[i] for i in np.argmax(
            alignment_score_query, axis=1)]
    alignment_score_query['pre_label'] = alignment_score_query['pre_label'].str.replace(
        '{}_'.format(sq_ref_SAMap), '').str.replace('{}_'.format(sq_que_SAMap), '')
    alignment_score_query['max_prob'] = alignment_score_query.filter(
        regex="^{}".format(sq_ref_SAMap)).max(axis=1)
    alignment_score_query.head(2)

    # constrect the res_obs
    res_obs = samap.adata.obs.loc[:, ['species']].copy()
    res_obs = res_obs.join(pd.DataFrame(samap.adata.obsm['X_umap'],
                                        columns='UMAP1,UMAP2'.split(','),
                                        index=samap.adata.obs.index))
    res_obs['dataset'] = res_obs['species'].map(
        {k: v for k, v in zip([sq_ref_SAMap, sq_que_SAMap], dsnames)})
    res_obs['cell_type'] = samap.adata.obs['{};{}_mapping_scores'.format(key_class1, key_class2)].str.replace(
        '^{}_'.format(sq_ref_SAMap), '', regex=True
    ).str.replace(
        '^{}_'.format(sq_que_SAMap), '', regex=True
    )
    res_obs['true_label'] = res_obs['cell_type']
    res_obs = res_obs.drop(columns='species')

    res_obs = res_obs.join(
        alignment_score_query.loc[:, 'pre_label,max_prob'.split(',')])
    res_obs['is_right'] = res_obs.eval('true_label == pre_label')
    res_obs.to_csv(resdir.joinpath('obs.csv'), index=True)

    # df_ratio
    df_ratio = group_agg(res_obs, 'dataset,is_right'.split(','), {
        'is_right': ['sum']
    }).merge(
        res_obs['dataset'].value_counts().to_frame(name='dataset_count'),
        on='dataset'
    )
    df_ratio = df_ratio.query('is_right').drop(
        columns='is_right').rename(
        columns={
            'dataset_count': 'is_right_count'})
    df_ratio = df_ratio.join(df_ratio['dataset'].str.extract(
        '(?P<tissue>[^_]+)_(?P<sp>[^_]+)'))
    assert df_ratio['sp'].is_unique, '[Error] not unique'
    df_ratio.index = df_ratio['sp'].to_numpy()
    df_ratio.loc[sp1, 'is_right_sum'] = np.nan

    df_ratio['ratio'] = df_ratio.eval("is_right_sum/is_right_count")

    df_ratio['type'] = 'species'
    df_ratio['name'] = ''
    df_ratio = df_ratio.loc[:,
                            'tissue,type,sp,name,is_right_sum,is_right_count,ratio'.split(',')]
    df_ratio.to_csv(resdir.joinpath('ratio.csv'), index=False)

    # plot umap
    samap.adata.obs['cell_type'] = samap.adata.obs['{};{}_mapping_scores'.format(key_class1, key_class2)].str.replace(
        '^{}_'.format(sq_ref_SAMap), '', regex=True).str.replace('^{}_'.format(sq_que_SAMap), '', regex=True)
    samap.adata.obs['dataset'] = samap.adata.obs['species'].map(
        {k: v for k, v in zip([sq_ref_SAMap, sq_que_SAMap], dsnames)})
    display(samap.adata.obs.head(2))

    ax = sc.pl.umap(samap.adata, color='dataset', show=False)
    ax.figure.savefig(
        Path(figdir).joinpath('umap_dataset.png'),
        bbox_inches="tight",
        dpi=120)
    ax = sc.pl.umap(samap.adata, color='cell_type', show=False)
    ax.figure.savefig(
        Path(figdir).joinpath('umap_umap.png'),
        bbox_inches="tight",
        dpi=120)

    # 完成标记
    resdir.joinpath("kvargs.json").write_text(dumps(kvargs))
    finish_content.append("[end] %f" % time.time())
    p_finish.write_text("\n".join(finish_content))


# ## run_cross_species_models

# In[ ]:


def run_cross_species_models(
    path_adata1,
    path_adata2,
    key_class1,
    key_class2,
    sp1,
    sp2,
    tissue_name,
    aligned=False,
    resdir_tag="",
    resdir=Path(
        "/public/workspace/licanchengup/download/test/test_result"
    ),
    limite_func=lambda adata1, adata2: (adata1, adata2),
    models='came,csMAHN'.split(','),
    **kvargs
):
    """
        kvargs:
        n_epochs:
            default,[100, 200, 300]
            stages,即res_0,res_1，res_2 的 epochs
            累加制，res_0,res_1，res_2,实际epochs分别为100,300,600
            故最终epochs为stages之和
            stages = kvargs.setdefault("n_epochs",[100, 200, 300])


        is_1v1: bool
            default,False
"""
    map_func = {
        'came': run_came,
        'csMAHN': run_csMAHN,
        'SAMap': run_SAMap
    }
    assert pd.Series([model in map_func.keys() for model in models]).all(
    ), "[Error] not all models in {}\nmodels {}".format(','.join(map_func.keys()), ','.join(models), )
    for model in models:
        path_varmap = get_path_varmap(
            map_sp[sp1], map_sp[sp2], model=model)
        # print('[path_varmap] {}\t{}'.format(model, Path(path_varmap).name))
        map_func[model](
            path_adata1,
            path_adata2,
            key_class1,
            key_class2,
            sp1,
            sp2,
            tissue_name,
            path_varmap,
            aligned=aligned,
            resdir_tag=";".join([model, resdir_tag]),
            resdir=resdir,
            limite_func=limite_func,
            **kvargs,
        )

