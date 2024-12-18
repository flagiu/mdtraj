#!/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 'large'
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import joblib

# NPT data
sc_args={"label":"SC", "color":"blue"}
a7_args={"label":"A7", "color":"red"}
T_eV=0.0345
P_eVang3=0.00624

ls=np.arange(4,8+1) #np.array([4,5,6,7,8],dtype=int)
ave=False  # use locally averaged q ?
plot=False
#

def plot_q_cross(folder, ls, axes, ave=True, nbins=30, color="black", label=None, contour_spacing="log"):
    print("Plotting q_l for many l's ...")
    nl=len(ls)
    for i1,l1 in enumerate(ls):
        ax=axes[i1,i1]
        ax.set_ylabel(r"$P(q)$")
        if ave:
            ax.set_xlabel(r"$\overline{q}_{%d}$"%l1)
            timestep,particle,q1 = np.loadtxt(folder+"/boo_ave.l%d.dat"%l1, unpack=True)
        else:
            ax.set_xlabel(r"$q_{%d}$"%l1)
            timestep,particle,q1 = np.loadtxt(folder+"/boo.l%d.dat"%l1, unpack=True)
        N=int(np.max(particle)+1)
        q1avg = q1.reshape(N,-1).mean(axis=0)
        h,edges=np.histogram(q1, bins=nbins, density=True)
        centers=(edges[1:]+edges[:-1])/2
        ax.plot(centers,h,alpha=0.9,color=color,linewidth=1)
        ax.axvline(q1avg.mean(), color=color, linewidth=2, linestyle="-")
        for i2 in range(i1+1,len(ls)):
            l2=ls[i2]
            ax=axes[i1,i2]
            if ave:
                ax.set_xlabel(r"$\overline{q}_{%d}$"%l1)
                ax.set_ylabel(r"$\overline{q}_{%d}$"%l2)
                timestep,particle,q2 = np.loadtxt(folder+"/boo_ave.l%d.dat"%l2, unpack=True)
            else:
                ax.set_xlabel(r"$q_{%d}$"%l1)
                ax.set_ylabel(r"$q_{%d}$"%l2)
                timestep,particle,q2 = np.loadtxt(folder+"/boo.l%d.dat"%l2, unpack=True)
            q2avg = q2.reshape(N,-1).mean(axis=0)
            h,xbins,ybins=np.histogram2d(q1,q2, bins=nbins, density=True)
            n=int(np.log10(h.max()))
            if contour_spacing=="log":
                contour_levels=np.array([0.1,0.7,0.9])*h.max()
            elif contour_spacing=="lin":
                contour_levels=np.linspace(h.min()*1.1,h.max()*0.9,3)
            else:
                print("Error: invalid keyword contour_spacing")
                sys.exit(1)
            X,Y = np.meshgrid(xbins[:-1],ybins[:-1])
            cs = ax.contour(X,Y, h.T, linewidths=0.5, levels=contour_levels, colors=color)
            ax.clabel(cs, fontsize=8)
            ax.scatter(q1avg,q2avg,marker='.',alpha=0.9,color=color,label=label)
    return
#-----------------------------------------------------------------------------------------#

if plot:
    fig, axes = plt.subplots(len(ls),len(ls), figsize=(30,30), dpi=300) #30,30
    fol="generate_sc/run_T%s_P%s/analysis"%(T_eV,P_eVang3)
    plot_q_cross(fol, ls, axes, ave=ave, **sc_args)
    fol="generate_a7/run_Seed389_T%s_P%s/analysis"%(T_eV,P_eVang3)
    plot_q_cross(fol, ls, axes, ave=ave, **a7_args)
    axes[0,0].legend(frameon=False, ncol=2)
    for ax in axes.ravel():
        ax.tick_params(which="both", direction="in")
    fig.tight_layout()
    fig.savefig("classification_cross_1GPa_A7_vs_SC.png")
    print("Saved: classification_cross_1GPa_A7_vs_SC.png")

#-------------------------------------------------------------------------------#

def classify(X, y, class_dict, axes, axes_predict, nbins=30, cm=plt.get_cmap("RdBu"), rseed=42 ,ax_confusionmatrix=None, plot=False):
    """
                 X: (input) nsamples x nfeatures
                 y: (ordered classes) nsamples
        class_dict: dictionary from 0,...,nclasses to text labels
              axes: nfeatures x nfeatures axes grid
      axes_predict: nfeatures x nfeatures axes grid
    """
    assert len(X.shape)==2
    assert len(y.shape)==1
    nsamples,nfeatures = X.shape
    assert len(y)==nsamples

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=rseed
    )

    def get_nice_hist(h,edges):
        bin=edges[1]-edges[0]
        centers=(edges[1:]+edges[:-1])/2
        np.insert(centers,0,edges[0]-bin/2)
        np.append(centers, edges[-1]+bin/2)
        hh=h.copy()
        np.insert(hh,0,0.)
        np.append(hh,0.)
        return hh,centers

    # Plot the training points
    if plot:
        edgecol_scatter = "black"
        for i in range(nfeatures):
            ax=axes[i,i]
            for cl in np.unique(y_train):
                sel=(y_train==cl)
                h,edges=np.histogram(X_train[sel,i],bins=nbins, density=True)
                h,centers=get_nice_hist(h,edges)
                ax.plot(centers,h,'-',alpha=0.9,color=cm(cl),linewidth=1,label="n.%d class %d"%(sum(sel),cl))
            for j in range(i+1,nfeatures):
                ax=axes[i,j]
                ax.scatter(X_train[:, i], X_train[:, j], c=y_train, marker='o',
                    cmap=cm, alpha=0.6, edgecolors=edgecol_scatter)



    clf = LogisticRegression(random_state=rseed, max_iter=10000)
    #clf = MLPClassifier(hidden_layer_sizes=(0,),solver='lbfgs', alpha=1e-5,
    #    max_iter=1000, random_state=rseed)

    clf = make_pipeline(StandardScaler(), clf)
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    # Plot the prediction on testing points
    if plot:
        sel_correct=(y_pred==y_test)
        for i in range(nfeatures):
            ax=axes_predict[i,i]
            for label in class_dict.keys():
                cl=class_dict[label]
                sel=(y_test==cl)&(sel_correct)
                h,edges=np.histogram(X_test[sel,i],bins=nbins, density=True)
                h,centers=get_nice_hist(h,edges)
                ax.plot(centers,h,'-',alpha=0.9,color=cm(cl),linewidth=1,label="%s correct"%(label))
                sel=(y_test==cl)&(~sel_correct)
                h,edges=np.histogram(X_test[sel,i],bins=nbins, density=True)
                h,centers=get_nice_hist(h,edges)
                ax.plot(centers,h,'--',alpha=0.9,color=cm(cl),linewidth=1,label="%s wrong"%(label))
            for j in range(i+1,nfeatures):
                ax=axes_predict[i,j]
                ax.scatter(
                    X_test[sel_correct, i], X_test[sel_correct, j], c=y_pred[sel_correct], marker='o', label="%correct",
                    cmap=cm, alpha=0.6, edgecolors=edgecol_scatter)
                ax.scatter(
                    X_test[~sel_correct, i], X_test[~sel_correct, j], c=y_pred[~sel_correct], marker='d', label="%wrong",
                    cmap=cm, alpha=0.6, edgecolors=edgecol_scatter)

    #cm = confusion_matrix(y_test, predictions, labels=clf.classes_)
    #cmdisp = ConfusionMatrixDisplay(confusion_matrix=cm,
    #                              display_labels=clf.classes_)
    SCtoSC=sum((y_test==class_dict["SC"])&(y_pred==y_test))
    SCtoA7=sum((y_test==class_dict["SC"])&(y_pred!=y_test))
    SC=SCtoSC+SCtoA7
    A7toA7=sum((y_test==class_dict["A7"])&(y_pred==y_test))
    A7toSC=sum((y_test==class_dict["A7"])&(y_pred!=y_test))
    A7=A7toA7+A7toSC
    print("SC-->SC:",SCtoSC)
    print("SC-->A7:",SCtoA7)
    print("A7-->A7:",A7toA7)
    print("A7-->SC:",A7toSC)
    confusion_summary=""
    confusion_summary+="         Accuracy: %.2f%%\n"%(100*(SCtoSC+A7toA7)/(SC+A7))
    confusion_summary+="SC, A7 guess rate: %.2f%%, %.2f%%\n"%(100*SCtoSC/SC, 100*A7toA7/A7)
    confusion_summary+="Balanced Accuracy: %.2f%%\n"%(100* (SCtoSC/SC + A7toA7/A7)/2 )
    confusion_summary+="SC, A7  miss rate: %.2f%%, %.2f%%\n"%(100*SCtoA7/SC, 100*A7toSC/A7)
    confusion_summary+="Balanced  missing: %.2f%%\n"%(100* (SCtoA7/SC + A7toSC/A7)/2 )
    print(confusion_summary)

    if plot:
        cmdisp = ConfusionMatrixDisplay.from_estimator(
            clf, X_test, y_test)
        cmdisp.plot(ax=ax_confusionmatrix)
        if ax_confusionmatrix is not None:
            ax_confusionmatrix.xaxis.set_ticklabels(class_dict.keys());
            ax_confusionmatrix.yaxis.set_ticklabels(class_dict.keys());
            ax_confusionmatrix.set_title(confusion_summary)
        else:
            plt.show()

    return clf


class_dict={"SC":0,"A7":1}
print(" Labels:",end="")
for cl in class_dict.keys():
    print("%s=%d "%(cl,class_dict[cl]),end="")
print()

Xa7 = []
for l in ls:
    if ave:
        timestep,particle,q = np.loadtxt("generate_a7/run_Seed389_T%s_P%s/analysis/boo_ave.l%d.dat"%(T_eV,P_eVang3,l), unpack=True)
    else:
        timestep,particle,q = np.loadtxt("generate_a7/run_Seed389_T%s_P%s/analysis/boo.l%d.dat"%(T_eV,P_eVang3,l), unpack=True)
    Xa7.append(q)
# shape: (N_a7*len_unique_timesteps_a7,len(ls))
Xa7=np.stack(Xa7,axis=1)
print(" A7 data shape:",Xa7.shape)


Xsc = []
for l in ls:
    if ave:
        timestep,particle,q = np.loadtxt("generate_sc/run_T%s_P%s/analysis/boo_ave.l%d.dat"%(T_eV,P_eVang3,l), unpack=True)
    else:
        timestep,particle,q = np.loadtxt("generate_sc/run_T%s_P%s/analysis/boo.l%d.dat"%(T_eV,P_eVang3,l), unpack=True)
    Xsc.append(q)
# shape: (N_sc*len_unique_timesteps_sc,len(ls))
Xsc=np.stack(Xsc,axis=1)
print(" SC data shape:",Xsc.shape," --> reduced to same as A7")
Xsc=Xsc[:Xa7.shape[0]]


Xdata = np.concatenate([Xsc,Xa7], axis=0)
Xlabel = np.empty(Xdata.shape[0])
Xlabel[:Xsc.shape[0]] = class_dict["SC"]
Xlabel[-Xa7.shape[0]:] = class_dict["A7"]

print(" Xdata shape:",Xdata.shape)
print(" Xlabel shape:",Xlabel.shape)

if plot:
    fig,axes = plt.subplots(len(ls),len(ls),figsize=(30,30),dpi=300)
    fig.suptitle("Train")

    fig_pred,axes_pred = plt.subplots(len(ls),len(ls),figsize=(30,30),dpi=300)
    fig_pred.suptitle("Test")

    fig_cm,ax_cm = plt.subplots(figsize=(8,8),dpi=300)
else:
    axes,axes_pred,ax_cm = None,None,None

clf = classify(Xdata,Xlabel, class_dict, axes,axes_pred,
    ax_confusionmatrix=ax_cm, plot=plot)
# save
joblib.dump(clf, "classifier.pkl")
# load
#clf2 = joblib.load("classifier.pkl")

if plot:
    fig_cm.tight_layout()
    fig_cm.savefig("classification_cm_1GPa_A7_vs_SC.png")
    print("Saved: classification_cm_1GPa_A7_vs_SC.png")

    for axx in [axes,axes_pred]:
        for i1 in range(len(ls)):
            ax=axx[i1,i1]
            l1=ls[i1]
            if ave:
                ax.set_xlabel(r"$\overline{q}_%d$"%l1)
                ax.set_ylabel(r"$P(\overline{q}_%d)$"%l1)
            else:
                ax.set_xlabel(r"${q}_%d$"%l1)
                ax.set_ylabel(r"$P({q}_%d)$"%l1)
            ax.legend(loc="center left", fontsize=8)
            ax.tick_params(which="both", direction="in")
            for i2 in range(i1+1,len(ls)):
                ax=axx[i1,i2]
                l2=ls[i2]
                if ave:
                    ax.set_xlabel(r"$\overline{q}_%d$"%l1)
                    ax.set_ylabel(r"$\overline{q}_%d$"%l2)
                else:
                    ax.set_xlabel(r"${q}_%d$"%l1)
                    ax.set_ylabel(r"${q}_%d$"%l2)
                ax.legend(loc="center left", fontsize=8)
                ax.tick_params(which="both", direction="in")

    fig.tight_layout()
    fig.savefig("classification_train_1GPa_A7_vs_SC.png")
    print("Saved: classification_train_1GPa_A7_vs_SC.png")

    fig_pred.tight_layout()
    fig_pred.savefig("classification_test_1GPa_A7_vs_SC.png")
    print("Saved: classification_test_1GPa_A7_vs_SC.png")

#-----------------------------------------------------#
