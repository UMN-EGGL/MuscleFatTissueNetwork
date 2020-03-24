import os
import glob
import re

import numpy as np
import pandas as pd
import camoco as co
import minus80 as m80
import pybiomart as pb
from bioservices.kegg import KEGG

from urllib.parse import urlparse
from collections import defaultdict,Counter

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt 


camoco_base = '/home/rob/.camoco/databases'

nets = ['EcAdipose','EcMuscle']

rule all:
    input:
        # Network Densities for KEGG
        expand("Data/KEGG/{net}_densities.csv",net=nets)
        # Network JSON for Cytoscape
        #expand('Health/{net}/{net}.json',net=nets)

        # Network Health Info (GO)
        #expand('Health/{net}/{net}_Health.summary.txt',net=nets),
        #expand('Health/{net}/{net}_Health_GO.csv',net=nets)

rule create_network_plot:
    input:
        EcMuscle = f"{camoco_base}/Expr.EcMuscle.db",
        EcAdipose = f"{camoco_base}/Expr.EcAdipose.db",
    output:
        fig = "Manuscript/Figures/Figure_1_Network_plots.eps"
    run:
        EcMuscle  = co.COB('EcMuscle')
        EcAdipose = co.COB('EcAdipose')

        # do degree distribution
        def degree_dist(cob,ax=None,fontsize=8):
            import powerlaw
            degree = cob.degree["Degree"].values
            # Using powerlaw makes run-time warning the first time you use it.
            # This is still an open issue on the creators github.
            # The creator recommends removing this warning as long as there is a fit.
            np.seterr(divide="ignore", invalid="ignore")
            fit = powerlaw.Fit(degree, discrete=True, xmin=1)
            # get an axis
            if ax is None:
                ax = plt.subplot()
            # Plot!
            emp = fit.plot_ccdf(ax=ax, color="k", linewidth=0.75, label="Empirical Data")
            pwr = fit.power_law.plot_ccdf(
                ax=ax, linewidth=1, color="k", linestyle=":", label="Power Law"
            )
            tpw = fit.truncated_power_law.plot_ccdf(
                ax=ax, linewidth=1, color="k", linestyle="-.", label="Truncated Power"
            )
            ax.set_ylabel("p(Degree≥x)", size=fontsize)
            ax.set_xlabel("Degree Frequency", size=fontsize)
           #plt.setp(ax.get_xticklabels(),size=fontsize)
           #plt.setp(ax.get_yticklabels(),size=fontsize)
            legend = ax.legend(
                loc="best",
                fontsize=fontsize-2,
                facecolor='white',
                edgecolor='white'
            )
            return ax

        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(7,4.5)) 
        node_size=10
        label_size=8
        cluster_line_width=1.5

        EcMuscle.plot_network(
            #target_genes=EcMuscle.MCL['MCL0'].loci,
            draw_clusters=True,min_cluster_size=20,cluster_std=3,
            node_size=node_size,min_degree=5,cluster_line_width=cluster_line_width,
            draw_edges=False,edge_alpha=0.1,max_edges=None,
            max_clusters=10,color_clusters=False,
            background_color='#888888',label_size=label_size,
            ax=ax1
        ) 
        EcAdipose.plot_network(
            #target_genes=EcMuscle.MCL['MCL0'].loci,
            draw_clusters=True,min_cluster_size=20,cluster_std=2,
            node_size=node_size,min_degree=5,cluster_line_width=cluster_line_width,
            draw_edges=False,edge_alpha=0.1,max_edges=None,
            max_clusters=10,color_clusters=False,
            background_color='#888888',label_size=label_size,
            ax=ax2
        ) 
        # Adjust the margins for the inset
        ax1.margins(0.01)
        ax2.margins(0.01)

        def adjust_margins(
                ax,
                top_scale=0,
                bottom_scale=0,
                right_scale=0,
                left_scale=0,
            ):
            # top

            # bottom
            bottom,top = ax.get_ylim()
            bottom_padding = abs(bottom - top) * bottom_scale
            top_padding = abs(bottom-top) * top_scale
            bottom = bottom - bottom_padding
            top = top + top_padding
            ax.set_ylim(bottom,top)
            # right
            left,right = ax.get_xlim()
            right_padding = abs(left - right) * right_scale
            left_padding = abs(left-right) * left_scale
            left = left - left_padding
            right = right + right_padding
            ax.set_xlim(left,right)

        adjust_margins(
            ax1,
            bottom_scale = 0.2,
            )

        adjust_margins(
            ax2,
            bottom_scale = 0,
            right_scale  = 0.0
        )

        # Insets
        EcMuscleInset = ax1.inset_axes([0,0,0.5,0.25],facecolor='white')
        degree_dist(EcMuscle,ax=EcMuscleInset)
        EcAdiposeInset = ax2.inset_axes([0,0,0.5,0.25],facecolor='white')
        degree_dist(EcAdipose,ax=EcAdiposeInset)

        plt.savefig(output.fig)

    
rule calculate_MCL_ontology_enrichment:
    input:
        EcGO = f'{camoco_base}/GOnt.EcGO.db',
        EcKEGG = f'{camoco_base}/Ontology.EcKEGG.db'
    output:
        table = 'Manuscript/Tables/Supplemental/Supplemental_Table_7_Network_MCL_Ontology_enrichment.xls'
    run:
        EcGO = co.GOnt('EcGO')
        EcKEGG = co.Ontology('EcKEGG')
        EcMuscle  = co.COB('EcMuscle')
        EcAdipose = co.COB('EcAdipose')
        # Create GO enrichments
        EcMuscleMCLGO = EcMuscle.MCL.enrichment(
            EcGO,
            min_term_size=10,
            max_term_size=10e10,
            return_table=True
        ).query(
            'bonferroni == True and num_common > 2 and target_term_size <= 300'
        ).sort_values(
            'num_common',
            ascending=False
        ).copy() 

        EcAdiposeMCLGO = EcAdipose.MCL.enrichment(
            EcGO,
            min_term_size=10,
            max_term_size=10e10,
            return_table=True
        ).query(
            'bonferroni == True and num_common > 2 and target_term_size <= 300'
        ).sort_values(
            'num_common',
            ascending=False
        ).copy() 

        # Create KEGG enrichments
        EcMuscleMCLKEGG = EcMuscle.MCL.enrichment(
            EcKEGG,
            min_term_size=10,
            max_term_size=10e10,
            return_table=True
        ).query(
            'bonferroni == True and num_common > 2 and target_term_size <= 300'
        ).sort_values(
            'num_common',
            ascending=False
        ).copy() 

        EcAdiposeMCLKEGG = EcAdipose.MCL.enrichment(
            EcKEGG,
            min_term_size=10,
            max_term_size=10e10,
            return_table=True
        ).query(
            'bonferroni == True and num_common > 2 and target_term_size <= 300'
        ).sort_values(
            'num_common',
            ascending=False
        ).copy() 

        # smash it all together
        MCL_Ont = pd.concat([
            EcMuscleMCLGO,
            EcMuscleMCLKEGG,   
            EcAdiposeMCLGO, 
            EcAdiposeMCLKEGG   
        ])

        MCL_Ont['Ontology'] = [x.split('_')[0] for x in MCL_Ont['label']]
        MCL_Ont['Term'] = [x.split('_')[1] for x in MCL_Ont['label']]

        MCL_Ont = MCL_Ont[[
            'source','name','source_term_size', 'Ontology',
            'Term','target_term_size','num_common',
            'num_universe','pval','bonferroni']
        ].copy()

        MCL_Ont.columns = [
             'Network','Cluster','Cluster_size','Ontology',
             'Term','Term_size','num_genes_common',
             'num_genes_universe','pval','bonferroni'
        ]
        MCL_Ont.to_excel(output.table)

        # These are the values for Table 2
        # ------------------------------------------- 
       #EcAdiposeMCLTerms = set([x.name for x in EcAdipose.MCL.terms(min_term_size=10)]) 
       #EcMuscleMCLTerms  = set([x.name for x in EcMuscle.MCL.terms(min_term_size=10)])
       #
       ## Calculate the overlaps
       #muscle_no_coclus  = set([x for x in EcMuscleMCLTerms if x not in mcl_enrich.source_term.unique()])
       #muscle_coclus     = set([x for x in EcMuscleMCLTerms if x in mcl_enrich.source_term.unique()])
       #muscle_no_ont     = EcMuscleMCLTerms.difference(MCL_Ont.query('Network == "EcMuscleMCL"').Cluster.unique())
       #muscle_ont        = EcMuscleMCLTerms.intersection(MCL_Ont.query('Network == "EcMuscleMCL"').Cluster.unique())

       #muscle_no_coclus_no_ont = muscle_no_coclus & muscle_no_ont
       #muscle_coclus_no_ont    = muscle_coclus & muscle_no_ont
       #muscle_no_coclus_ont    = muscle_no_coclus & muscle_ont
       #muscle_coclus_ont       = muscle_coclus & muscle_ont

       #adipose_no_coclus  = set([x for x in EcAdiposeMCLTerms if x not in mcl_enrich.target_term.unique()])
       #adipose_coclus     = set([x for x in EcAdiposeMCLTerms if x in mcl_enrich.target_term.unique()])
       #adipose_no_ont     = EcAdiposeMCLTerms.difference(MCL_Ont.query('Network == "EcAdiposeMCL"').Cluster.unique())
       #adipose_ont        = EcAdiposeMCLTerms.intersection(MCL_Ont.query('Network == "EcAdiposeMCL"').Cluster.unique())

       #adipose_no_coclus_no_ont = adipose_no_coclus & adipose_no_ont
       #adipose_coclus_no_ont    = adipose_coclus & adipose_no_ont
       #adipose_no_coclus_ont    = adipose_no_coclus & adipose_ont
       #adipose_coclus_ont       = adipose_coclus & adipose_ont
    
rule calculate_cross_MCL_enrichment:
    input:
        EcMuscle = f"{camoco_base}/Expr.EcMuscle.db",
        EcAdipose = f"{camoco_base}/Expr.EcAdipose.db"
    output:
        table = 'Manuscript/Tables/Supplemental/Supplemental_Table_6_Network_MCL_co-clusters.xls'
    run:
        EcMuscle  = co.COB('EcMuscle')
        EcAdipose = co.COB('EcAdipose')
        
        # Create the table --------------------------------------------
        # Do enrichment
        mcl_enrich = EcMuscle.MCL.enrichment(
            EcAdipose.MCL,
            min_term_size=10,
            max_term_size=10e10,
            return_table=True
        ) 
        mcl_enrich['target_network'] = [x.split('_')[0] for x in mcl_enrich['label']]
        mcl_enrich['target_term'] = [x.split('_')[1] for x in mcl_enrich['label']]
        mcl_enrich = mcl_enrich.query('bonferroni == True and num_common > 2')

        mcl_enrich = mcl_enrich[
            ['source','name','source_term_size','target_network',
             'target_term','target_term_size','num_common',
             'num_universe','pval','bonferroni']
        ]
        mcl_enrich.columns = [
             'source_network','source_term','source_term_size','target_network',
             'target_term','target_term_size','num_common',
             'num_universe','pval','bonferroni'
        ]
        mcl_enrich.to_excel(output.table,index=False)
        

rule plot_mcl_upset:
    input:
        table = 'Manuscript/Tables/Supplemental/Supplemental_Table_6_Network_MCL_co-clusters.xls'
    output:
        fig = "Manuscript/Figures/Figure_3_EcMuscle_EcAdipose_MCL_co-clusters.eps" 
    run:
        # Read in the table
        mcl_enrich = pd.read_excel(input.table)
        # generate some data structures
        EcMuscle  = co.COB('EcMuscle')
        EcAdipose = co.COB('EcAdipose')

        EcAdiposeMCLTerms = set([x.name for x in EcAdipose.MCL.terms(min_term_size=10)]) 
        EcAdiposeMCLTerms = sorted(EcAdiposeMCLTerms,key=lambda x: int(x.replace('MCL','')))
        EcMuscleMCLTerms = set([x.name for x in EcMuscle.MCL.terms(min_term_size=10)])
        EcMuscleMCLTerms = sorted(EcMuscleMCLTerms,key=lambda x: int(x.replace('MCL','')),reverse=True)

        def upset_plot(
            x,
            y,
            intersect,
            xlabel,
            ylabel
        ):
            # Upset plot
            fig = plt.figure(facecolor='black',figsize=(7,9))
            gs = fig.add_gridspec(
                14,16,
                hspace=0.1,
                wspace=0.1
            ) 
            hax = fig.add_subplot(gs[2:,2:])
            ax1 = fig.add_subplot(gs[:2,2:],sharex=hax)
            ax2 = fig.add_subplot(gs[2:,:2],sharey=hax)

            xindex = np.arange(0,len(x),1)
            yindex = np.arange(0,len(y),1)
            xlabels = {x:i for i,x in enumerate(x)}
            ylabels = {y:i for i,y in enumerate(y)}

            xcoor = []
            ycoor = []

            for x,y in intersect:
                xcoor.append(xlabels[x])
                ycoor.append(ylabels[y])

            hax.scatter(
                xcoor,
                ycoor,
                color='k',
                marker='D'
            )        
            hax.set_xticks(xindex)
            hax.set_xticklabels(
                labels=xlabels.keys(),
                rotation=45,
                ha='right',
                rotation_mode="anchor",
                fontsize='xx-small', 
            )
            
            for label in hax.get_xticklabels()[-2::-2]:
                label.set_visible(False)
            hax.set_yticks(np.arange(0,len(ylabels)))
            hax.set_xlabel(xlabel) 

            hax.yaxis.tick_right()
            hax.set_yticklabels(
                labels=ylabels.keys(),
                fontsize='xx-small'
            )
            for label in hax.get_yticklabels()[1::2]:
                label.set_visible(False)

            hax.set_ylabel(ylabel)
            hax.yaxis.set_label_position("right") 
            
            # EcAdipose bar plot
            x_coclusters = defaultdict(int,Counter([x[0] for x in intersect]))
            ax1.set_title("Network MCL co-clusters")
            ax1.bar(
                x=xindex,
                height=[x_coclusters[x] for x in xlabels.keys()],
                color='k'
            )
            ax1.grid(False,axis='x')
            plt.setp(ax1.get_xticklabels(),visible=False)
            ax1.set_yticks(
                np.arange(5,max(x_coclusters.values()),5),
            )
            ax1.set_ylabel('Num. \nco-clusters')
            # EcMuscle bar plot
            y_coclusters = defaultdict(int,Counter([x[1] for x in intersect]))
            x = np.arange(len(y_coclusters))
            ax2.barh(
                y=yindex,
                width=[-1*y_coclusters[y] for y in ylabels],
                color='k'
            )
            ax2.yaxis.set_ticks_position("right") 
            ax2.grid(False,axis='y')
            ax2.set_xticks(
                np.arange(-2,-1*max(y_coclusters.values()),-2)
            )
            ax2.set_xticklabels(
                labels=np.arange(2,max(y_coclusters.values()),2)
            )
            ax2.set_xlabel('Num. \nco-clusters')

            # Final cleanup
            hax.margins(0.01)
            plt.subplots_adjust(
                left=0.02
            )
            #plt.setp(ax2.get_yticklines(),visible=False)
            plt.setp(ax2.get_yticklabels(),visible=False)
            return fig

        upset_plot(
            EcMuscleMCLTerms,
            EcAdiposeMCLTerms,
            list(zip(mcl_enrich.source_term,mcl_enrich.target_term)),
            'EcMuscle MCL',
            'EcAdipose MCL'
        )

        # Save figure
        plt.savefig(output.fig)


rule calculate_kegg_enrichment:
    input:
        f"{camoco_base}/Ontology.EcKEGG.db"
    output:
        dens_out = f"Data/KEGG/{{net_name}}_densities.csv"
    params:
        kegg_name = "EcKEGG",
        num_bs = 1000
    run:        
        KEGG = co.Ontology(params.kegg_name)
        cob = co.COB(wildcards.net_name)
        df = pd.DataFrame(columns=["KeggTerm","desc","size","density","density_pval"])
        for term in KEGG:
            # this is necessary because we want to know the number of genes 
            # IN THE NETWORK, not the term
            term.loci = list(filter(lambda x: x in cob, term.loci))
            # SHOULD THIS BE term.loci ???
            if len(term) < 10 or len(term) > 300:
                continue
            # Calculate the emprirical score
            density = cob.density(term.loci)
            # Calculate the bootstrap score using the built-in refgen included
            # in the cob object this will guarantee that the random genes are 
            # actually in the network (as opposed to getting a random set from
            # EquCab3)
            bs_density = [cob.density(cob.refgen.random_genes(n=len(term.loci))) for _ in range(params.num_bs)]
            # calculate proportion of bs densities greater than the empirical density
            pval = sum(i > density for i in bs_density)/len(bs_density)
            # append the results to a data frame for use later
            df.loc[len(df)] = [term.name, 
                               term.desc.replace(" ","_"),
                               len(term.loci), 
                               density, 
                               pval]
            # A COUPLE TERM DESCRIPTIONS ARE QUOTED WHICH IS ANNOYING ME...
            df.to_csv(output.dens_out, index = False)


rule compare_network_GO_terms:
    input:
        muscle_go = 'Health/EcMuscle/EcMuscle_Health_GO.csv',
        adipose_go = 'Health/EcAdipose/EcAdipose_Health_GO.csv',
        muscle_kegg = 'Data/KEGG/EcMuscle_densities.csv',
        adipose_kegg = 'Data/KEGG/EcAdipose_densities.csv',
        go=f'{camoco_base}/GOnt.EcGO.db'
    output:
        go_xls  = "Manuscript/Tables/Supplemental/Supplemental_Table_4_GOTerm_Densities.xlsx",
        kegg_xls= "Manuscript/Tables/Supplemental/Supplemental_Table_5_KEGGTerm_Densities.xlsx",
        figure  = "Manuscript/Figures/Figure_2_GO_Term_Network_comparison.eps"
    run:
        import camoco as co
        import pandas as pd
        import numpy as np
        import matplotlib
        import PIL # not needed but must be installed to save to tiff
        from matplotlib import pylab as plt
        # Build data structures
        EcGO = co.GOnt('EcGO')
        mgo = pd.read_table(input.muscle_go,sep=',')
        mkegg = pd.read_table(input.muscle_kegg,sep=',',dtype={'KeggTerm':str})
        fgo = pd.read_table(input.adipose_go,sep=',')
        fkegg = pd.read_table(input.adipose_kegg,sep=',',dtype={'KeggTerm':str})
        # Filter out negative densities
        mgo.density.values[mgo.density < 0] = 0
        mgo.density_pval.values[mgo.density < 0] = 1
        fgo.density.values[fgo.density < 0] = 0
        fgo.density_pval.values[fgo.density < 0] = 1
        # same for KEGG 
        mkegg.density.values[mkegg.density < 0] = 0
        mkegg.density_pval.values[mkegg.density < 0] = 1
        fkegg.density.values[fkegg.density < 0] = 0
        fkegg.density_pval.values[fkegg.density < 0] = 1
        # Merge KEGG terms
        kegg = mkegg[['KeggTerm','desc','size','density','density_pval']].merge(
            fkegg[['KeggTerm','density','density_pval']],
            how='inner',
            on='KeggTerm'
        ).set_index('KeggTerm')
        kegg.rename(
            columns={
                'density_x' : 'muscle_density',
                'density_pval_x':'muscle_density_pval',
                'density_y' : 'adipose_density',
                'density_pval_y':'adipose_density_pval',
            },
            inplace=True
        )
        # Merge the p-values on terms that have data for both "inner"
        go = mgo[["GOTerm",'desc','size','density','density_pval']].merge(
            fgo[["GOTerm",'density',"density_pval"]],
            how='inner',
            on="GOTerm"
        ).set_index('GOTerm')
        go.rename(
            columns={
                'density_x' : 'muscle_density',
                'density_pval_x':'muscle_density_pval',
                'density_y' : 'adipose_density',
                'density_pval_y':'adipose_density_pval',
            },
            inplace=True
        )       
        # add in a columns with a namespace mapping
        namespace_map = {k:v for k,v in EcGO.db.cursor().execute('SELECT term,val FROM term_attrs WHERE key = "namespace"')}
        go.insert(0,'namespace',[namespace_map.get(x,'missing') for x in go.index.values])
        # add in the name column because we forgot it :( (use SQL because its WAY faster)
        name_map = {k:v for k,v in  EcGO.db.cursor().execute("SELECT id,name FROM terms").fetchall()}
        go.insert(1,'name',[name_map[x] for x in go.index.values])

        # Calculate ben-hoch FDR -------------------------------------------------------------------------
        FDR=0.01
        go.sort_values('muscle_density_pval',inplace=True)
        go['muscle_bh'] = [p <= (k/len(go))*FDR  for k,p in enumerate(go.muscle_density_pval)]
        go.sort_values('adipose_density_pval',inplace=True)
        go['adipose_bh'] = [p <= (k/len(go))*FDR  for k,p in enumerate(go.adipose_density_pval)]

        kegg.sort_values('muscle_density_pval',inplace=True)
        kegg['muscle_bh'] = [p <= (k/len(kegg))*FDR  for k,p in enumerate(kegg.muscle_density_pval)]
        kegg.sort_values('adipose_density_pval',inplace=True)
        kegg['adipose_bh'] = [p <= (k/len(kegg))*FDR  for k,p in enumerate(kegg.adipose_density_pval)]

        # convert to -log10 pval
        go['muscle_density_minus_log_pval'] = -1*np.log10(go.muscle_density_pval + 0.001)
        go['adipose_density_minus_log_pval'] = -1*np.log10(go.adipose_density_pval + 0.001)

        kegg['muscle_density_minus_log_pval'] = -1*np.log10(kegg.muscle_density_pval + 0.001)
        kegg['adipose_density_minus_log_pval'] = -1*np.log10(kegg.adipose_density_pval + 0.001)
   
        # Figure out the -log10 pvalue for FDR 1%
        go_bh_FDR = np.mean([
            min(go.muscle_density_minus_log_pval[go.muscle_bh]),
            min(go.adipose_density_minus_log_pval[go.adipose_bh])
        ])
        kegg_bh_FDR = np.mean([
            min(kegg.muscle_density_minus_log_pval[kegg.muscle_bh]),
            min(kegg.adipose_density_minus_log_pval[kegg.adipose_bh])
        ])
         
        SMALL_SIZE   = 8
        MEDIUM_SIZE  = 10
        BIGGER_SIZE  = 12
        BIGGEST_SIZE = 14
        
        fig,((g1,k1),(g2,k2),(g3,k3)) = plt.subplots(
            3,2,
            figsize=(7,7),
            #constrained_layout=True
            tight_layout=True,
            gridspec_kw={
                'hspace': 0.55,
                'wspace': 0.22
            }
        )
        plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
       #plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
       #plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
       #plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
       #plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels

        plt.rc('axes', titlesize=BIGGEST_SIZE)    

        plt.margins(0.0)
        #plt.tight_layout()
        # First Panel ----------------------------------------------------
        # EcMuscle score vs pval
        x = g1.scatter(
            go.muscle_density,
            go.muscle_density_minus_log_pval,
            color='k',
            s=10
        )
        l = g1.hlines(
            go_bh_FDR,
            min(go.muscle_density),
            max(go.muscle_density),
            color='darkgray'
        )
        g1.set_xlabel('EcMuscle Density Score')
        g1.set_yticks([0,1,2,3])
        g1.set_ylabel('-log10(p-value)')
        g1.legend(['GO Term'],loc='lower right')
        g1.set_title('GO')
        # Do KEGG
        x = k1.scatter(
            kegg.muscle_density,
            kegg.muscle_density_minus_log_pval,
            color='k',
            s=10
        )
        l = k1.hlines(
            kegg_bh_FDR,
            min(kegg.muscle_density),
            max(kegg.muscle_density),
            color='darkgray'
        )
        k1.set_xlabel('EcMuscle Density Score')
        k1.set_yticks([0,1,2,3])
        k1.legend(['KEGG Term'],loc='lower right')
        k1.set_title('KEGG')
                                                                                        
        # Second Panel ----------------------------------------------------
        # EcAdipose score vs pval 
        g2.scatter(
            go.adipose_density,
            go.adipose_density_minus_log_pval,
            color='k',
            s=10
        )
        g2.hlines(
            go_bh_FDR,
            min(go.adipose_density),
            max(go.adipose_density),
            color='darkgray'
        )
        g2.set_xlabel('EcAdipose Density Score')
        g2.set_yticks([0,1,2,3])
        g2.set_ylabel('-log10(p-value)')
        g2.legend(['GO Term'],loc='lower right')
        # Do KEGG
        k2.scatter(
            kegg.adipose_density,
            kegg.adipose_density_minus_log_pval,
            color='k',
            s=10
        )
        k2.hlines(
            kegg_bh_FDR,
            min(kegg.adipose_density),
            max(kegg.adipose_density),
            color='darkgray'
        )
        k2.annotate(
            'FDR 1%',
            xy=(max(kegg.muscle_density)-50,kegg_bh_FDR+5),
            xycoords='data'
        )
        k2.set_xlabel('EcAdipose Density Score')
        k2.set_yticks([0,1,2,3])
        k2.legend(['KEGG Term'],loc='lower right')

        # Third Panel ----------------------------------------------------
        # Get numbers for the shared and specific values
        both_go = [
            len(go.query(f'muscle_density_minus_log_pval >= {go_bh_FDR} and adipose_density_minus_log_pval >= {go_bh_FDR}')),
            len(go.query(f'muscle_density_minus_log_pval >= {go_bh_FDR} and adipose_density_minus_log_pval >= {go_bh_FDR}'))  
        ]
        tissue_only_go = [
            len(go.query(f'muscle_density_minus_log_pval >= {go_bh_FDR} and adipose_density_minus_log_pval < {go_bh_FDR}')),
            len(go.query(f'muscle_density_minus_log_pval < {go_bh_FDR} and adipose_density_minus_log_pval >= {go_bh_FDR}'))
        ]
        both_kegg = [
            len(kegg.query(f'muscle_density_minus_log_pval >= {kegg_bh_FDR} and adipose_density_minus_log_pval >= {kegg_bh_FDR}')),
            len(kegg.query(f'muscle_density_minus_log_pval >= {kegg_bh_FDR} and adipose_density_minus_log_pval >= {kegg_bh_FDR}'))  
        ]        
        tissue_only_kegg = [
            len(kegg.query(f'muscle_density_minus_log_pval >= {kegg_bh_FDR} and adipose_density_minus_log_pval < {kegg_bh_FDR}')),
            len(kegg.query(f'muscle_density_minus_log_pval < {kegg_bh_FDR} and adipose_density_minus_log_pval >= {kegg_bh_FDR}'))
        ]

        ind = [0,1]
        # plot the number of both first
        b1=g3.bar(
            ind,both_go,color='white',
            edgecolor='black',hatch='//',
            width=0.8
        )
        # then plot the next bar ontop of the previes
        b2=g3.bar(
            ind,tissue_only_go,bottom=both_go,
            color='white',edgecolor='black',hatch='\\\\',
            width=0.8
        )
        # set the y ticks to be the total GO Terms
        g3.set_xlim(-0.5,1.5)
        g3.set_xticks([0,1])
        g3.set_xticklabels(['EcMuscle','EcAdipose']) 
        g3.set_yticks(np.arange(0, 5418, 1000))
        g3.set_ylabel('Number Terms')
        #g3.legend((b1[0], b2[0]), ('Both Tissues', 'Single Tissue'))

        def go_to_percent(x):
            return 100*(x/len(go))

        def go_to_number(x):
            return (x*len(go))/100

        secg3 = g3.secondary_yaxis('right', functions=(go_to_percent, go_to_number))        
        secg3.set_yticks(np.arange(0,100,20))

        b1=k3.bar(
            ind,both_kegg,color='white',
            edgecolor='black',hatch='//', 
            width=0.8
        )
        # then plot the next bar ontop of the previes
        b2=k3.bar(
            ind,tissue_only_kegg,bottom=both_kegg,
            color='white',edgecolor='black',hatch='\\\\',
            width=0.8
        )

        # set the y ticks to be the total GO Terms
        k3.set_xlim(-0.5,1.5)
        k3.set_xticks([0,1])
        k3.set_xticklabels(['EcMuscle','EcAdipose']) 
        k3.set_yticks(np.arange(0, 70, 10))
        k3.set_yticklabels(['0','','20','','40','','60',''])
        k3.legend((b2[0], b1[0]), ('Single Tissue', 'Both Tissues'))

        def kegg_to_percent(x):
            return 100 * (x/len(kegg))

        def kegg_to_number(x):
            return (x * len(kegg))/100

        seck3 = k3.secondary_yaxis('right', functions=(kegg_to_percent, kegg_to_number))        
        seck3.set_yticks(np.arange(0,100,20))
        seck3.set_ylabel('% Terms')
 
        fig.savefig(output.figure)
        go.to_excel(output.go_xls)
        kegg.to_excel(output.kegg_xls)
        

rule generate_cytoscape_networks:
    input:
        network=f'{camoco_base}/Expr.{{net_name}}.db'
    output:
        json=f'Health/{{net_name}}/{{net_name}}.json'
    params:
        net_name=f'{{net_name}}'
    run:
        import camoco as co
        net = co.COB(params.net_name)
        net.to_json(
            filename=output.json,
            #max_edges=100000,
            ontology=net.MCL
        )

rule calculate_health:
    input:
        network=f'{camoco_base}/Expr.{{net_name}}.db'
    output:
        summary=f'Health/{{net_name}}/{{net_name}}_Health.summary.txt',
        go_scores=f'Health/{{net_name}}/{{net_name}}_Health_GO.csv'
    params:
        network=f'{{net_name}}'
    shell:
        '''
            camoco \
            health \
            {params.network} \
            --go EcGO \
            --refgen EquCab3 \
            --num-bootstraps 1000 \
            --out-dir Health/{params.network}/
        '''    

rule build_fat_net:
    input:
        fat='Data/Expr/Equine_Adipose_TPM.tsv',
        refgen=f'{camoco_base}/RefGen.EquCab3.db'
    output:
        f'{camoco_base}/Expr.EcAdipose.db'
    shell:
        '''
            camoco \
            build-cob \
            {input.fat} \
            EcAdipose \
            "EcAdipose Network" \
            EquCab3
        '''

rule build_muscle_net:
    input:
        muscle='Data/Expr/Equine_Muscle_TPM.tsv',
        refgen=f'{camoco_base}/RefGen.EquCab3.db'
    output:
        f'{camoco_base}/Expr.EcMuscle.db'
    shell:
        '''
            camoco \
            build-cob \
            {input.muscle} \
            EcMuscle \
            "EcMuscle Network" \
            EquCab3 \
        '''

rule build_kegg:
    input:
        refgen = f"{camoco_base}/RefGen.EquCab3.db",
        kegg_gene_map = "Data/KEGG/EquCab3_KEGG.map.tsv"
    output:
        f"{camoco_base}/Ontology.EcKEGG.db"
    run:
        EquCab3 = co.RefGen(input.refgen)
        kegg_d = defaultdict(list)
        with open(input.kegg_gene_map) as f:
            next(f)
            for line in f:
                line = line.strip().split("\t")
                # get kegg pathway id and gene id 
                pw = line[1].split("+")[0]
                gene_id = line[-1]
                # dict with pathway:RefGen gene locus object
                kegg_d[pw].append(EquCab3.from_id(gene_id))
        kegg_terms = []
        # generate list of terms with kegg names and pathway descriptions
        for k, v in kegg_d.items():
            pw_desc = KEGG().find("pathway",k).strip().split("\t")[1]
            if pw_desc:
                print(f"Found KEGG pathway for {k} ({pw_desc})")
            else:
                print(f"Unable to find pathway for {k}")
            # THOUGHTS ABOUT TO BETTER TO DEAL WITH THIS POTENTIAL SITUATION???
            kegg_terms.append(co.Term(k, desc=pw_desc, loci=v))
        # kegg ontology from kegg terms above    
        KEGG = co.Ontology.from_terms(name="EcKEGG",
                                      description="EquCab3 KEGG ontology",
                                      refgen=EquCab3,
                                      terms=kegg_terms) 

rule build_go:
    input:
        gene_mappings = 'Data/GO/EquCab3_GO.tsv',
        obo = 'Data/GO/go.obo',
        refgen=f'{camoco_base}/RefGen.EquCab3.db'
    output:
        f'{camoco_base}/GOnt.EcGO.db'
    shell:
        '''
            camoco \
            build-go \
            {input.gene_mappings} \
            {input.obo} \
            "EcGO" \
            "EquCab3 Gene Ontology" \
            EquCab3 \
            --go-col 1 \
            --id-col 3 
        ''' 

rule build_go2:
    input:
        gene_mappings = 'Data/GO/EquCab3_GO.map.tsv',
        obo = 'Data/GO/go.obo',
        refgen=f'{camoco_base}/RefGen.EquCab3.db'
    output:
        f'{camoco_base}/GOnt.EcGO2.db'
    shell:
        '''
            camoco \
            build-go \
            {input.gene_mappings} \
            {input.obo} \
            "EcGO2" \
            "EquCab3 Gene Ontology 2" \
            EquCab3 \
            --go-col 1 \
            --id-col 3 
        ''' 

rule download_go_obo:
    output:
        obo_file = 'Data/GO/go.obo'
    shell:
        '''
            curl http://current.geneontology.org/ontology/go.obo -o {output.obo_file}
        '''

rule generate_ontology_mappings:
    input:
        f'{camoco_base}/RefGen.EquCab3.db',
    output:
        go_out = "Data/GO/EquCab3_GO.map.tsv",
        kegg_out = "Data/KEGG/EquCab3_KEGG.map.tsv",
        unknown = "RefGen_Unknown.map.tsv"
    run:
        # Create an ID map from ncbi xref gene IDS to the gene1, gene2, gene3 ...
        id_map = {}
        for g in co.RefGen("EquCab3").iter_genes():
            # create dict first split on comma, then split on first colon
            xrefs = {k: v for k, v in [x.split(":", 1) for x in g.attr["Dbxref"].split(",")]}
            # now xrefs is a k:v mapping between db name and ID (e.g. {'GeneID': '100067156', 'VGNC': 'VGNC:23791'})
            id_map[xrefs["GeneID"]] = g.attr["ID"]

        # Get ensembl biomart dataset and query ncbi, go, and kegg attributes
        dataset = pb.Dataset(
            name="ecaballus_gene_ensembl", 
            host="http://www.ensembl.org"
        )
        ecab_df = dataset.query(
            attributes=[
                "ensembl_gene_id", "go_id",
                "entrezgene_id",   "kegg_enzyme"
            ]
        )
        ecab_df.columns = ["ensecag", "go", "ncbi", "kegg"]
        # convert ecab dataframe to dictionary
        dedup_df = ecab_df.drop_duplicates(
            subset=["ensecag", "go", "ncbi", "kegg"]
        ).fillna("").to_dict("index")
        
        in_refgen = {}
        not_refgen = {}
        # sort biomart entries by presence/absence in refgen db
        for k, v in dedup_df.items():
            if v["ncbi"]:
                v["ncbi"] = str(int(v["ncbi"]))
                ncbi_attr = v["ncbi"]
                if ncbi_attr in id_map:
                    in_refgen[k] = v
                    in_refgen[k]["geneid"] = id_map[ncbi_attr]
                else:
                    not_refgen[k] = v
            else:
                not_refgen[k] = v
        # write all entries to file based on go, kegg, or missing from refgen db
        with open(output.go_out, "w") as go_out, \
             open(output.kegg_out, "w") as kegg_out, \
             open(output.unknown, "w") as unknown :
            print("ensecag", "go", "ncbi", "geneid", sep="\t", file=go_out)
            print("ensecag", "kegg", "ncbi", "geneid", sep="\t", file=kegg_out)
            print("ensecag", "go", "kegg", "ncbi", sep="\t", file=unknown)
            for k, v in in_refgen.items():
                if v["go"]:
                    print(v["ensecag"], 
                          v["go"],
                          v["ncbi"],
                          v["geneid"], 
                          sep="\t", file=go_out)
                if v["kegg"]:
                    print(v["ensecag"],
                          v["kegg"], 
                          v["ncbi"], 
                          v["geneid"], 
                          sep="\t", file=kegg_out)
    
            for k, v in not_refgen.items():
                print(v["ensecag"], 
                      v["go"], 
                      v["kegg"], 
                      str(v["ncbi"]), 
                      sep="\t", file=unknown)

 
rule build_refgen:
    input:
        gff='Data/RefGen/GCF_002863925.1_EquCab3.0_genomic.nice.gff.gz'
    output:
        f'{camoco_base}/RefGen.EquCab3.db'
    shell:
        '''
            camoco \
            build-refgen \
            {input.gff} \
            "EquCab3" \
            "EquCab3 nicificed" \
        '''
            
