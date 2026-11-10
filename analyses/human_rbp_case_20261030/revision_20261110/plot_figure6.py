"""Figure 6: cohort profiles, cross-cell correspondence and FAM120A follow-up.

Read-only reuse of the original cohort outputs and the audited FAM120A data.
The previous publication figure is never overwritten.
"""
from pathlib import Path
import sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.collections import PolyCollection
from matplotlib.ticker import FuncFormatter

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
OLD=ROOT/'revision_20261030'
FOLLOW=ROOT/'fam120a_followup_20261109'
OUT=HERE/'publication';OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(OLD));sys.path.insert(0,str(ROOT/'scripts'))
from plot_revision import profile_matrix, annotate_functions, paired_distances, FUNCTIONS, REGIONS, REGION_LABELS
from plot_integration import logo
from plot_common import COLORS,CELL,DENSITY

plt.rcParams.update({'font.family':'Arial','font.size':7.0,'axes.labelsize':7.0,
 'axes.titlesize':7.8,'xtick.labelsize':6.1,'ytick.labelsize':6.1,'pdf.fonttype':42,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,
 'legend.fontsize':6.1,'legend.frameon':False,'xtick.major.width':.5,'ytick.major.width':.5})
REGCOL={'UTR3':COLORS['lavender'],'Intron':COLORS['green'],'Other':'#E7E5EB'}
TRACKS=[('ENCFF318HXI','K562','IP 1'),('ENCFF161DFB','K562','IP 2'),
 ('ENCFF135HYO','K562','Input'),('ENCFF354DCE','HepG2','IP 1'),
 ('ENCFF354PCR','HepG2','IP 2'),('ENCFF341HEH','HepG2','Input')]

def label(fig,ax,letter,title,pad=8):
    box=ax.get_position()
    fig.text(box.x0-.035,box.y1+pad/fig.bbox.height*fig.dpi/72,letter,
             fontsize=10.5,fontweight='bold',va='bottom')
    ax.set_title(title,loc='left',pad=pad)

def save(fig,name):
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    overflow=[]
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text().strip():continue
        bb=text.get_window_extent(renderer)
        # Text attached to cropped axes can be outside the page only when hidden.
        if bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5:
            overflow.append(text.get_text())
    if overflow:raise ValueError('Text outside figure: '+repr(overflow))
    for ext in ['png','pdf']:fig.savefig(OUT/(name+'.'+ext),dpi=600,facecolor='white')
    plt.close(fig)

def main_figure():
    d=pd.read_csv(ROOT/'data/landscape/datasets.tsv',sep='\t').sort_values(
        ['Intron','protein','cell'],ascending=[False,True,True]).reset_index(drop=True)
    profiles=np.load(ROOT/'summary/profiles.npz')
    tests=pd.read_csv(OLD/'data/all_function_tests.tsv',sep='\t')
    pairs=pd.read_csv(OLD/'data/regional_pairs.tsv',sep='\t')
    overlap=pd.read_csv(OLD/'data/coordinate_overlap.tsv',sep='\t')
    report=json.loads((OLD/'data/report.json').read_text())
    d[['accession','protein','cell','n_filtered','Intron']].to_csv(OUT/'profile_row_order.tsv',sep='\t',index=False)
    fig=plt.figure(figsize=(183/25.4,245/25.4))
    a=fig.add_axes([.080,.943,.895,.030]);a.axis('off')
    label(fig,a,'A','Large-scale positional profiling',pad=7)
    boxes=[(0,.31,'223 eCLIP datasets','150 RBPs · 844854 IDR peaks',COLORS['pale_blue']),
           (.345,.30,'12 RNAmeta analyses','2676 completed module outputs',COLORS['lavender']),
           (.68,.32,'Two human cell lines','K562 and HepG2',COLORS['pink'])]
    for x,w,title,subtitle,color in boxes:
        a.add_patch(Rectangle((x,0),w,.90,facecolor=color,alpha=.42,lw=0,transform=a.transAxes))
        a.text(x+w/2,.59,title,ha='center',va='center',fontsize=7.7,fontweight='bold',color=COLORS['deep'])
        a.text(x+w/2,.23,subtitle,ha='center',va='center',fontsize=6.0)
    axes=[fig.add_axes([x,.703,.185,.193]) for x in [.095,.33,.56,.79]]
    b=axes[0];label(fig,b,'B','Region composition',pad=17)
    im=b.imshow(d[REGIONS],aspect='auto',cmap=DENSITY,vmin=0,vmax=1,
        extent=[-.5,5.5,223.5,.5],interpolation='nearest')
    b.set_xticks(range(6),REGION_LABELS,rotation=45,ha='right');b.tick_params(axis='x',length=0,pad=3)
    for edge in np.arange(.5,5.5):b.axvline(edge,color='white',lw=.5)
    b.set_yticks([1,50,100,150,200,223]);b.set_ylabel('Dataset rank',labelpad=2)
    b.tick_params(axis='y',length=2,pad=2)
    strip=b.inset_axes([1.025,0,.026,1]);strip.imshow(
        np.array([matplotlib.colors.to_rgba(CELL[c]) for c in d.cell])[:,None,:],
        aspect='auto',interpolation='nearest');strip.axis('off')
    b.legend(handles=[Patch(facecolor=CELL[c],label=c) for c in ['K562','HepG2']],
        loc='center',bbox_to_anchor=(.5,-.35),ncol=2,handlelength=.9,columnspacing=.7,
        borderpad=0,borderaxespad=0,fontsize=6)
    cmap=DENSITY.copy();cmap.set_bad('#D0D0D0')
    for ax,letter,key,title in zip(axes[1:],'CDE',
      ['metagene_mrna','splicesite__5PSS','splicesite__3PSS'],
      ['mRNA distribution','5′ splice boundary','3′ splice boundary']):
        label(fig,ax,letter,title,pad=17)
        matrix=profile_matrix(d,profiles,key);x=profiles[key+'__x']
        if key.startswith('metagene'):matrix=matrix*3
        limits=[0,3] if key.startswith('metagene') else [x[0],x[-1]]
        hh=ax.imshow(matrix,aspect='auto',cmap=cmap,vmin=0,vmax=4,
              extent=[*limits,223.5,.5],interpolation='nearest');ax.set_yticks([])
        if key=='metagene_mrna':
            ax.set_xticks([])
            for xx in [1,2]:ax.axvline(xx,color='white',lw=.5)
            for j,(name,color) in enumerate(zip(['5′UTR','CDS','3′UTR'],
                [COLORS['pale_blue'],COLORS['green'],COLORS['lavender']])):
                ax.add_patch(Rectangle((j,-.10),1,.075,transform=ax.get_xaxis_transform(),
                    facecolor=color,edgecolor='white',lw=.5,alpha=.6,clip_on=False))
                ax.text(j+.5,-.0625,name,transform=ax.get_xaxis_transform(),ha='center',va='center',fontsize=6)
            ax.set_xlabel('Normalized transcript position',fontsize=6,labelpad=19)
        else:
            ax.set_xticks(np.arange(np.ceil(x[0]/100)*100,x[-1]+1,100))
            ax.axvline(0,color=COLORS['deep'],ls='--',lw=.5);ax.set_xlabel('Distance (nt)',labelpad=2)
            names=['Exon','Intron'] if key.endswith('5PSS') else ['Intron','Exon']
            for lo,hi,name in zip([x[0],0],[0,x[-1]],names):
                ax.plot([lo,hi],[1.009]*2,transform=ax.get_xaxis_transform(),color=COLORS['deep'],lw=.5,clip_on=False)
                ax.text((lo+hi)/2,1.018,name,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=6)
        for spine in ax.spines.values():spine.set_visible(False)
        ax.spines['bottom'].set_visible(True)
    for spine in ['left','top','right']:b.spines[spine].set_visible(False)
    cb=fig.colorbar(im,cax=fig.add_axes([.095,.618,.185,.006]),orientation='horizontal')
    cb.set_label('Fraction of peaks',fontsize=6,labelpad=1);cb.ax.tick_params(labelsize=5.8,length=2,pad=1)
    cb=fig.colorbar(hh,cax=fig.add_axes([.33,.618,.645,.006]),orientation='horizontal')
    cb.set_label('Relative peak density',fontsize=6,labelpad=1);cb.ax.tick_params(labelsize=5.8,length=2,pad=1)
    f=fig.add_axes([.195,.448,.155,.113]);label(fig,f,'F','Functional associations')
    hm=annotate_functions(f,tests,FUNCTIONS,['UTR5','CDS','UTR3','Intron'])
    f.tick_params(axis='y',labelsize=5.9,pad=2);f.tick_params(axis='x',labelsize=5.9)
    cb=fig.colorbar(hm,cax=fig.add_axes([.364,.448,.007,.113]),ticks=[-1,0,1])
    cb.set_label('Cliff’s δ',fontsize=6,labelpad=1);cb.ax.tick_params(labelsize=5.8,length=2,pad=1)
    g=fig.add_axes([.465,.448,.205,.113])
    paired_distances(g,pairs,'Cross-cell correspondence','G',report['regional_cross_cell'])
    g.texts[0].remove();g.set_title('');label(fig,g,'G','Cross-cell correspondence')
    g.set_ylabel('Distributional distance',labelpad=2);g.tick_params(axis='x',labelsize=6)
    h=fig.add_axes([.780,.448,.195,.113]);label(fig,h,'H','Coordinates and distributions')
    for flag,color in [(True,COLORS['blue']),(False,COLORS['pink'])]:
        sub=overlap[(overlap.difference>0)==flag]
        h.scatter(sub.mean_peak_overlap*100,sub.same_distance,s=10,color=color,alpha=.8,edgecolors='white',linewidth=.25)
    fam=overlap[overlap.protein.eq('FAM120A')].iloc[0]
    h.scatter([fam.mean_peak_overlap*100],[fam.same_distance],s=27,facecolors='none',edgecolors=COLORS['deep'],linewidth=.8,zorder=4)
    h.annotate('FAM120A',(fam.mean_peak_overlap*100,fam.same_distance),xytext=(9,4),textcoords='offset points',
               fontsize=6.2,color=COLORS['deep'],arrowprops=dict(arrowstyle='-',color=COLORS['deep'],lw=.5))
    h.set_xlim(-2,102);h.set_ylim(0,.65);h.set_xticks([0,50,100]);h.set_xlabel('Mean peak overlap (%)',labelpad=2)
    h.set_ylabel('Same-RBP distance',labelpad=2)
    h.text(.04,.98,'69 / 73 retain closer profiles',va='top',transform=h.transAxes,fontsize=5.8)
    # New FAM120A evidence, displayed without changing any prior calculations.
    summary=pd.read_csv(FOLLOW/'data/regional_expression_sensitivity.tsv',sep='\t')
    i=fig.add_axes([.135,.262,.345,.117]);label(fig,i,'I','FAM120A regional composition')
    cohorts=[('all_peaks','All peaks'),('baseline_tpm1_fc2','Expression matched'),('stable_common_bound_min1','Shared targets')]
    for k,(cohort,name) in enumerate(cohorts):
        for j,cell in enumerate(['K562','HepG2']):
            sub=summary[(summary['mode']=='center_sense')&(summary.cohort==cohort)&(summary.cell==cell)]
            s=sub.set_index('region').fraction;vals=[s.UTR3,s.Intron,1-s.UTR3-s.Intron];left=0;y=k*2.75+j
            for color,value in zip(REGCOL.values(),vals):
                i.barh(y,value*100,left=left*100,height=.66,color=color,edgecolor='white',linewidth=.3)
                if value>.11:i.text((left+value/2)*100,y,f'{value*100:.1f}',ha='center',va='center',fontsize=5.7)
                left+=value
            i.text(103,y,str(int(sub.n_peaks.iloc[0])),va='center',fontsize=5.8)
        i.text(0,k*2.75-.77,name,fontsize=6.1,ha='left')
    i.set_ylim(7.1,-1.5);i.set_xlim(0,123);i.set_xticks([0,50,100]);i.set_xlabel('Peaks (%)',fontsize=6.2,labelpad=2)
    i.set_yticks([k*2.75+j for k in range(3) for j in range(2)],['K562','HepG2']*3)
    i.tick_params(axis='y',length=0,pad=3,labelsize=6);i.spines['left'].set_visible(False);i.spines['bottom'].set_bounds(0,100)
    i.text(107,-1.1,'n',fontsize=6,ha='center')
    i.legend(handles=[Patch(facecolor=c,label=name) for c,name in zip(REGCOL.values(),['3′UTR','Intron','Other regions'])],
        loc='upper left',bbox_to_anchor=(-.02,-.30),ncol=3,columnspacing=.7,handlelength=.9,fontsize=5.8,borderpad=0)
    jax=fig.add_axes([.675,.274,.300,.105]);label(fig,jax,'J','Within shared targets')
    paired=pd.read_csv(FOLLOW/'data/paired_gene_comparison.tsv',sep='\t')
    paired=paired[paired.min_peaks_per_cell.eq(1)]
    for yy,region in zip([1,0],['UTR3','Intron']):
        s=paired[paired.region.eq(region)].iloc[0];x=100*s.mean_difference_HepG2_minus_K562
        jax.errorbar(x,yy,xerr=[[x-100*s.bootstrap95_low],[100*s.bootstrap95_high-x]],fmt='o',ms=4.3,
           color=REGCOL[region],lw=1.4,capsize=2)
        jax.text(x,yy+.23,f'{x:+.1f}',fontsize=6.3,ha='center',color=COLORS['deep'])
    jax.axvline(0,color='#AAA5B2',ls=':',lw=.6);jax.set_xlim(-29,26);jax.set_ylim(-.6,1.75)
    jax.set_yticks([1,0],['3′UTR','Intron']);jax.set_xticks([-20,0,20]);jax.spines['left'].set_visible(False)
    jax.tick_params(axis='y',length=0);jax.set_xlabel('HepG2 − K562 (percentage points)',fontsize=6.2,labelpad=2)
    jax.text(.98,1.0,'n = 354 genes',transform=jax.transAxes,ha='right',va='bottom',fontsize=6.1)
    # Complete locus: four IP replicates and both inputs on the same scale.
    signal=pd.read_csv(FOLLOW/'data/lmnb2_signal.tsv.gz',sep='\t')
    peaks=pd.read_csv(FOLLOW/'data/annotated_peaks_with_expression.tsv.gz',sep='\t')
    peaks=peaks[peaks['mode'].eq('center_sense')&peaks.gene_id.eq('ENSG00000176619')]
    model=pd.read_csv(FOLLOW/'data/lmnb2_transcript_features.tsv',sep='\t')
    model=model[model.transcript.eq('ENST00000325327.3')]
    lo,hi=2427600,2457500; ymax=45
    axes=[]
    for row,(acc,cell,name) in enumerate(TRACKS):
        ax=fig.add_axes([.15,.170-row*.019,.825,.015]);axes.append(ax)
        sub=signal[signal.accession.eq(acc)]
        vertices=[[(max(start,lo),0),(max(start,lo),abs(v)),(min(end,hi),abs(v)),(min(end,hi),0)]
                  for start,end,v in zip(sub.start,sub.end,sub.score)]
        ax.add_collection(PolyCollection(vertices,facecolors=CELL[cell] if name!='Input' else '#AAA6B1',
                          edgecolors='none',rasterized=True))
        ax.set_xlim(lo,hi);ax.set_ylim(0,ymax);ax.set_xticks([]);ax.set_yticks([])
        for spine in ax.spines.values():spine.set_visible(False)
        ax.axhline(0,color='#DBD8E0',lw=.35)
        ax.text(-.015,.35,f'{cell} {name}',transform=ax.transAxes,ha='right',va='center',fontsize=5.8,
                color=CELL[cell] if name!='Input' else '#66616E')
    label(fig,axes[0],'K','FAM120A signal at LMNB2',pad=9)
    axes[0].text(1,1.1,'ENCODE signal: 0-45',transform=axes[0].transAxes,ha='right',va='bottom',fontsize=5.8)
    for row,cell in enumerate(['K562','HepG2']):
        ax=fig.add_axes([.15,.063-row*.009,.825,.006])
        for r in peaks[peaks.cell.eq(cell)].itertuples():
            ax.add_patch(Rectangle((r.start,0),r.end-r.start,1,facecolor=CELL[cell],lw=0))
        ax.set_xlim(lo,hi);ax.set_ylim(0,1);ax.axis('off')
        ax.text(-.015,.5,cell+' peaks',transform=ax.transAxes,ha='right',va='center',fontsize=5.8)
    ax=fig.add_axes([.15,.033,.825,.013]);exons=model[model.kind.eq('exon')]
    ax.hlines(0,exons.start.min(),exons.end.max(),color=COLORS['deep'],lw=.55)
    for r in exons.itertuples():ax.add_patch(Rectangle((r.start,-.16),r.end-r.start,.32,color=COLORS['deep'],lw=0))
    for r in model[model.kind.eq('CDS')].itertuples():ax.add_patch(Rectangle((r.start,-.3),r.end-r.start,.6,color=COLORS['deep'],lw=0))
    for pos in np.linspace(exons.start.min(),exons.end.max(),10)[1:-1]:ax.plot(pos,0,marker='<',ms=2,color=COLORS['deep'])
    ax.set_xlim(lo,hi);ax.set_ylim(-.5,.5);ax.set_yticks([]);ax.spines['left'].set_visible(False)
    ax.set_xticks([2428000,2436000,2444000,2452000]);ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x/1e6:.3f}'))
    ax.tick_params(axis='x',length=2,pad=1);ax.set_xlabel('chr19 (Mb), hg19',fontsize=6.1,labelpad=1)
    ax.text(-.015,0,'LMNB2',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=5.8)
    save(fig,'fig6')
    (OUT/'figure6_panel_map.json').write_text(json.dumps(dict(width_mm=183,height_mm=245,
      panels={'A':'cohort scale','B-E':'four aligned distribution views','F':'functional associations',
              'G':'cross-cell correspondence','H':'coordinate overlap and FAM120A identification',
              'I':'FAM120A regional composition','J':'354 paired genes, primary analysis',
              'K':'LMNB2 original eCLIP/input signal'},signal_scale=[0,45],
      signal_transform='absolute value of negative-strand score',
      gene_model='ENST00000325327.3',source_existing_figure=str(OLD/'publication/fig6.png')),indent=2))

def motif_supplement():
    d=pd.read_csv(ROOT/'data/landscape/datasets.tsv',sep='\t')
    profiles=np.load(ROOT/'summary/profiles.npz')
    motifs={m['accession']:m for m in json.loads((ROOT/'summary/motifs.json').read_text())}
    key='splicesite__3PSS';ids=profiles[key+'__ids'].tolist();x=profiles[key+'__x']
    fig=plt.figure(figsize=(183/25.4,112/25.4))
    for col,protein in enumerate(['QKI','PTBP1','U2AF2']):
        left=.09+col*.315;ax=fig.add_axes([left,.57,.215,.27]);label(fig,ax,'ABC'[col],protein)
        for cell in ['K562','HepG2']:
            rec=d[(d.protein==protein)&(d.cell==cell)].iloc[0]
            ax.plot(x,profiles[key][ids.index(rec.accession)],color=CELL[cell],lw=1.2,label=cell)
        ax.axvline(0,color=COLORS['deep'],lw=.6,ls='--');ax.set_xlim(-200,100);ax.set_xticks([-200,0,100])
        ax.set_xlabel('Distance to 3′ boundary (nt)',labelpad=2);ax.set_ylabel('Relative density',labelpad=2)
        ax.legend(loc='center',bbox_to_anchor=(.5,-.55),ncol=2,fontsize=6.2,handlelength=1.2,columnspacing=1)
        for row,cell in enumerate(['K562','HepG2']):
            rec=d[(d.protein==protein)&(d.cell==cell)].iloc[0];m=motifs[rec.accession]['motifs'][0]
            la=fig.add_axes([left,[.255,.055][row],.215,.115]);logo(la,m);width=len(m['pwm']['A'])
            la.set_xticks(np.arange(width)+.5,range(1,width+1));la.tick_params(axis='x',length=2,pad=1.5)
            la.spines['bottom'].set_visible(True);la.set_ylabel('Bits',labelpad=1)
            la.set_title(cell,loc='left',color=CELL[cell],fontsize=6.5,pad=2)
            if row==1:la.set_xlabel('Motif position',fontsize=6.5,labelpad=2)
    save(fig,'supplementary_fig17_rbp_boundary_motifs')

if __name__=='__main__':
    main_figure();motif_supplement()
    print('WROTE',OUT/'fig6.png')
