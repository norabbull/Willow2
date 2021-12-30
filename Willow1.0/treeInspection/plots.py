# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 06:45:58 2021

@author: norab
"""

from treeHelpers import *
import pandas as pd
from plotnine import *
import glob

#%%
# Load data 

# Folder with all files
folder = 'C:/Users/norab/Master/thesis_data/'

# Load data using helper functions
nz_phydists = load_nz_phydists(file=folder + 'meta_data/nonzero_phydists.csv')  # works
uniqseqs = load_uniqseqs(folder + 'meta_data/9381_uniqseqs.txt') #works
GDRs = load_GDRs(file = folder + 'result_data/GDR/GDR_all.csv') # works
GDRnull = load_GDRnull()   # works
GDRsuper = load_GDRs(file = folder + 'result_data/GDR/GDR_all.csv', group_category = 'super')
GDRsub = load_GDRs(file = folder + 'result_data/GDR/GDR_all.csv', group_category = 'sub')

# Merge data into common df
data_all = nz_phydists.merge(uniqseqs, on = 'gene')
data_all = data_all.merge(GDRs, on = 'gene')
data_all.drop_duplicates(inplace=True)

# Signifiacnt genes
data_all_sig_SUB = pd.read_csv('C:/Users/norab/Master/thesis_data/result_data/GDRsignificant/SDR1005_significantGenesSUB.csv')
data_all_sig_SUPER = pd.read_csv('C:/Users/norab/Master/thesis_data/result_data/GDRsignificant/SDR1005_significantGenesSUPER.csv')

# all_data, super
data_all_super = nz_phydists.merge(uniqseqs, on = 'gene')
data_all_super = data_all_super.merge(GDRsuper, on = 'gene')
data_all_super.drop_duplicates(inplace=True)

# all_data, sub
data_all_sub = nz_phydists.merge(uniqseqs, on = 'gene')
data_all_sub = data_all.merge(GDRsub, on = 'gene')
data_all.drop_duplicates(inplace=True)


#%% Stats

# Get basic statistics of GDR values
GDRs.describe()
GDRsuper.describe()
GDRsub.describe()


#%%

# =============================================================================
# Plot: GDR distributions(figure 6.5 b)
# Violin + interconnecting lines between GDR for super and sub
# =============================================================================



shift = 0.1

def alt_sign(x):
    "Alternate +1/-1 if x is even/odd"
    return (-1) ** x

m1 = aes(x=stage('level', after_scale='x+shift*alt_sign(x)'))              # shift outward
m2 = aes(x=stage('level', after_scale='x-shift*alt_sign(x)'), group='gene')  # shift inward



(ggplot(GDRs, aes('level', 'GDR', fill = 'level'))
 + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.8, show_legend = False)
 + geom_point(m2, color='none', alpha=0.6, size=1.5, show_legend=False)
 + geom_line(m2, color='saddlebrown', size=0.4, alpha=0.3)
 + scale_fill_manual(values=['indigo', 'darkorange'])
 + theme_classic()  # change color
 + theme(figure_size=(8, 6))
 + labs(title='GDR distributions')
 + xlab("Group level")
 + ylab("GDR")
 + theme_classic()
 + theme(
     axis_text_x=element_text(rotation=30, hjust=1),
         plot_title=element_text(color = "black",
                             size = 35,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
          axis_title=element_text(color='black',
                             size=25,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
          axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=25,
     legend_key_height=15,
     legend_key_size=25,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
     )
)



# =============================================================================
# Plot: GDR distributions (figure 6.5 a) 
# density + lines between sub and super
# =============================================================================


(ggplot(GDRs, aes('level', 'GDR', fill = 'level'))
 + geom_violin(m1, style = 'left-right', alpha = 0.7, size = 0.9, show_legend = False)
 + geom_boxplot(width = shift, alpha=0.4, size = 0.9, show_legend = False)
 + scale_fill_manual(values=['indigo', 'darkorange'])
 + theme_classic()  # change color
 + theme(figure_size=(8, 6))
 + labs(title='GDR distributions')
 + xlab("Group level")
 + ylab("GDR")
 + theme(
     axis_text_x=element_text(rotation=30, hjust=1),
         plot_title=element_text(color = "black",
                             size = 35,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
          axis_title=element_text(color='black',
                             size=25,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
          axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=25,
     legend_key_height=15,
     legend_key_size=25,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
     )
)



# =============================================================================
# Plot: GDR distributions (figure 6.4)
# Density plot 
# =============================================================================


(ggplot(GDRs, aes(x='GDR', fill='level'))
 + geom_density(adjust = 1/2, alpha=0.5) 
 + scale_fill_manual(values=['indigo', 'darkorange'])
 + theme_classic()  # change color
 + theme(figure_size=(8, 6))
 + labs(title='GDR distributions')
 + xlab("GDR")
 + ylab("Density")
 + theme(
     axis_text_x=element_text(rotation=30, hjust=1),
         plot_title=element_text(color = "black",
                             size = 35,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
          axis_title=element_text(color='black',
                             size=25,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
          axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=30,
     legend_key_height=20,
     legend_key_size=40,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold', size=12),
     legend_position='right'

     )
)



#%% Null distributions

# =============================================================================
# # GDR null distribution for 30.000 GDR null values (figure 6.6)
# =============================================================================

l = [pd.read_csv(filename) for filename in glob.glob("C:/Users/norab/Master/thesis_data/result_data/GDRnull/super/*.csv")]
all_GDRnull = pd.concat(l, axis=0)
all_GDRnull.describe()

half_GDRnull = all_GDRnull.iloc[0:30000]
half_GDRnull.columns = ['GDRnull']


(
ggplot(half_GDRnull, aes(x='GDRnull', y=after_stat('density')))
+ geom_histogram(
    colour='darkgreen', # change the outline
    size=2,        # change the thickness of the outline
    alpha=0.7      # change the transparency
    )
+ labs(title='30.000 random GDR values (super)')
+ scale_x_continuous(name="GDR")
+ scale_y_continuous(name="density")
  + theme(
     plot_title=element_text(color = "black",
                             size = 32,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
     axis_title=element_text(color='black',
                             size=22,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
     
     axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=30,
     legend_key_height=15,
     legend_key_size=30,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
     #legend_title_align='center',
)
  )



#%%

# =============================================================================
# Non - zero values in distance matrices (figure 6.8)
# =============================================================================

nz_phydists = nz_phydists.dropna()

(
ggplot(nz_phydists, aes(x='totdist', y=after_stat('count/sum(count)*100')))
+ geom_histogram(
    colour='darkgreen', # change the outline
    size=2,        # change the thickness of the outline
    alpha=0.7      # change the transparency
    )
+ labs(title='Non - zero values in distance matrices')
 + scale_x_continuous(name="% non-zero distance-values in distance matrices")
 + scale_y_continuous(name="% of all distance matrices")
  + theme(
     plot_title=element_text(color = "black",
                             size = 32,
                             family = 'serif',
                             weight = 'semibold',
                             margin={'b':20}),
     axis_title=element_text(color='black',
                             size=22,
                             family='sans-serif',
                             weight='bold',
              
                              margin={'t':15, 'r':15}),
     
     axis_text=element_text(size=12, weight='semibold'),
     legend_key_width=30,
     legend_key_height=15,
     legend_key_size=30,
     legend_entry_spacing=10,
     legend_box_margin=5,
     legend_title=element_text(weight='bold')
     #legend_title_align='center',
)
  )
