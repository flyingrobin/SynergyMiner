__author__ = 'liwang'

import pandas as pd
import numpy as np
import pymysql
import re
import networkx as nx
from collections import defaultdict


# get drug target list
def drugTargets_dict(df_drug_list):
    # get a dictionary in which the key is the drug name and the value is a list of target genes
    target_dict=dict()
    for index, row in df_drug_list.iterrows():
        targets=[re.sub('[\*\s]','',i) for i in row['Target_annotated'].split(',')]
        target_dict[row['ChallengeName']]=targets
    return target_dict

def get_target_combID(in_network_list,dict_target):

    # df: the data frame of drug combination data (containing synergy scores)
    # return a dictionary where key is the drug combination ID and the key is a list of all target combination IDs
    target_combIDs=defaultdict()

    for key_i in dict_target.keys():
        for key_j in dict_target.keys():
            if key_i in in_network_list and key_j in in_network_list and key_i != key_j:
                drug_combID=key_i+'.'+key_j
                list1=dict_target[key_i]
                list2=dict_target[key_j]

                target_combIDs[drug_combID]=[]
                for i in list1:
                    for j in list2:
                        target_combIDs[drug_combID].append([re.sub('[\s\*]','',i),re.sub('[\s\*]','',j)])

    return target_combIDs

# Extract features from the graph

def in_same_cluster(cluster_list,node2):
    for cluster in cluster_list:
        if node2 in cluster:
            return 1

    return 0

def neighbor_overlap_orderK(G,node1,node2,k):
    nei1=nx.single_source_shortest_path(G,node1,cutoff=k).keys()
    nei2=nx.single_source_shortest_path(G,node2,cutoff=k).keys()

    # return (x1, x2)
    # where x1 = the number of overlapped neighboring nodes
    # and   x2 = the number of overlapped neighboring nodes / the number of union(two neighborhood nodes)

    shared_neighbor=set(nei1).intersection(set(nei2))
    unioned_neighbor=set(nei1).union(set(nei2))

    return (len(shared_neighbor),float(len(shared_neighbor))/float(len(unioned_neighbor)))

def iter_feature1(G,node1,node2):
    # Extract graph features for a given pair of nodes
    # return a data frame of features


    feat_out=[]
    try:
        # Length of the shortest path
        feat = nx.shortest_path_length(G,node1,node2)
        feat_out.append(feat)

        # Number of shortest path
        feat = len(list(nx.all_shortest_paths(G,node1,node2)))
        feat_out.append(feat)

    except nx.NetworkXNoPath:
        feat_out.extend([0,0])

    # first order neighbood overlap
    feat1,feat2 = neighbor_overlap_orderK(G,node1,node2,1)
    feat_out.extend([feat1,feat2])

    # second order neighbood overlap
    feat1,feat2 = neighbor_overlap_orderK(G,node1,node2,2)
    feat_out.extend([feat1,feat2])

    # average neighbor degree
    #feat1,feat2 = nx.average_neighbor_degree(G,nodes=[node1,node2]).values()
    #feat_out.extend([feat1,feat2])

    # Connectivity
    feat = nx.node_connectivity(G,node1,node2)
    feat_out.append(feat)

    # whether the nodes are in the same cluster
    feat = in_same_cluster(nx.cliques_containing_node(G,node1),node2)
    feat_out.append(feat)

    return feat_out

def graph_features1(G,targetID_list):
    # return data frame
    feat_list=[]
    for node1,node2 in targetID_list:
        if node1 in G and node2 in G:
            feat_list.append(iter_feature1(G,node1,node2))

    df_feat=np.apply_along_axis(np.mean,0,feat_list)

    return df_feat

def key_iter_feat2(target_dict,dic):
    value_list=[]
    value2_list=[]
    index_list=[]

    df_feat=np.array([999,998])
    for comID in target_dict.keys():
        for target_pair in target_dict[comID]:
            node1,node2 = target_pair
            if node1 in dic and node2 in dic :
                value_list.append([dic[node1],dic[node2]])

        index_list.append(comID)
        df_feat = np.vstack((df_feat,np.apply_along_axis(np.mean,0,value_list)))

    return  pd.DataFrame(df_feat,index=['Header']+index_list)

def graph_features2(G, target_dict):

    # Run the graph calculation only once
    # Degree
    dict_degree = G.degree()

    # Centrality
    # degree centrality
    dict_degree_cen = nx.degree_centrality(G)

    # closeness centrality
    dict_close_cen = nx.closeness_centrality(G)

    # betweeness centrality
    dict_betwn_cen = nx.betweenness_centrality(G)


    # Clustering
    # clustering coefficient for the nodes
    dict_clst_coef = nx.clustering(G_targets_in_net)

    df_feat=pd.DataFrame()

    # Start iterating


    for each_dict in [dict_degree,dict_degree_cen,dict_close_cen,dict_betwn_cen,dict_clst_coef]:
        df_feat = pd.concat([df_feat,key_iter_feat2(target_dict,each_dict)],axis=1)

    return df_feat

def makePPIdata(in_network_list,df_drug_list,G):
    # Convert drug combination IDs to target gene IDs
    target_dict = drugTargets_dict(df_drug_list)
    target_combIDs=get_target_combID(in_network_list,target_dict)

    x_feat2=graph_features2(G,target_combIDs)
    print 'Done with feature set 2 extraction!'

    x_feat1=np.empty((0,8))
    idxID=[]
    count = 0
    for key in target_combIDs.keys():
        feat1=graph_features1(G,target_combIDs[key])
        print key,x_feat1.shape,feat1.shape
        x_feat1=np.vstack((x_feat1,[feat1]))
        idxID.append(key)
        count += 1
        if count % 50 == 0 :
            print 'Finished parsing',count,'combinations of feature set 1.'

    return (x_feat1,x_feat2)

def tabFromSQL(mydb,table):
        query_colname="SELECT Column_name from INFORMATION_SCHEMA.COLUMNS where TABLE_NAME='"+table+"';"
        query_tab="SELECT * FROM Insight."+table+";"
        with mydb:
            cur = mydb.cursor()
            cur.execute(query_colname)
            colnames=[i[0] for i in cur.fetchall()]

            cur.execute(query_tab)
            query_results = cur.fetchall()
        return pd.DataFrame(list(query_results),columns=colnames)

if __name__=="__main__":
    mydb = pymysql.connect(host='localhost',user='root',password='',db='Insight')

    # Read PPI data from MySQL database
    df_PPI_all=tabFromSQL(mydb,'PPI_all')
    df_PPI_map=tabFromSQL(mydb,'PPI_IDmap')

    # Create mapping dictionary, there are 41 Ids that map to multiple gene names so we can ignore them
    # For multi-mapping IDs, we take the last entry
    uniprot2gene=dict()
    for index, row in df_PPI_map.iterrows():
        uniprot2gene[row['From']]= row['To']

    # Read KEGG table from MySQL database
    df_KEGG=tabFromSQL(mydb,'KEGG_pathway')

    # Combine all kepp pathway nodes
    kegg_nodes=[item for sub in [df_KEGG[col] for col in df_KEGG.columns] for item in sub]
    kegg_nodes=list(set(kegg_nodes))

    # read Drug info data from MySQL database
    df_drug_list=tabFromSQL(mydb,'Drug_info')

    # Build PPI network graph
    G=nx.Graph()
    for index, row in df_PPI_all.iterrows():
        if row.ProteinA in uniprot2gene and row.ProteinB in uniprot2gene:
            G.add_edge(uniprot2gene[row.ProteinA],uniprot2gene[row.ProteinB])

    # generate graph containing both cancer pathway genes and target genes (union of both)
    target_dict = drugTargets_dict(df_drug_list)
    target_list = [item for sublist in target_dict.values() for item in sublist]
    target_in_net=list(set(kegg_nodes+target_list))
    G_targets_in_net=G.subgraph(target_in_net)

    # Check drugs that has no target genes in the network, essentially they are chemotherapy drugs
    in_network=[]
    for key,val in target_dict.items():
        not_network = True
        for gene in val:
            gene=re.sub('[\s\*]','',gene)
            if gene in G_targets_in_net:
                not_network = False
        if not not_network:
            in_network.append(key)

    x_feat1_2,x_feat2=makePPIdata(in_network,df_drug_list,G_targets_in_net)
    x_combID=pd.concat([x_feat1,x_feat2.drop(x_feat2.index[0])],axis=1)