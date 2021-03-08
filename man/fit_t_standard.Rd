\name{fit_t_standard}
\alias{fit_t_standard}

\title{
Fits standard models of trait evolution incorporating known and nuisance measurement error
}
\description{
Fits Brownian motion (BM), Ornstein-Uhlenbeck (OU), or early burst (EB) models of trait evolution to a given dataset and phylogeny.
}
\usage{
fit_t_standard(phylo, data, model=c("BM","OU","EB"), error, two.regime=FALSE, 
		method="Nelder-Mead", echo=TRUE, ...)
}

\arguments{
  \item{phylo}{
an object of type 'phylo' (see ape documentation); if \code{two.regime=TRUE}, this must be a simmap object from \code{\link{make.simmap}} with two regimes
}
  \item{data}{
a named vector of trait values with names matching \code{phylo$tip.label}
}
  \item{error}{
if including measurement error, a named vector of error values (the standard error of the [intraspecific] mean) with names matching \code{phylo$tip.label}
}

  \item{model}{
model chosen to fit trait data, \code{"BM"} is the Brownian motion model, \code{"OU"} is the Ornstein-Uhlenbeck model, and \code{"EB"} is the early burst model.
}
  \item{pars}{
vector specifying starting parameter values for maximum likelihood optimization. If unspecified, default values are used (see Details)
}

  \item{echo}{
prints information to console during fit
}

}

\value{
	a list with the following elements:
  \item{LH}{maximum log-likelihood value}
  \item{aic}{Akaike Information Criterion value}
  \item{aicc}{AIC value corrected for small sample size}
  \item{free.parameters}{number of free parameters from the model}
  \item{sig2}{maximum-likelihood estimate of \code{sig2} parameter}
  \item{alpha}{maximum-likelihood estimate of \code{alpha} parameter of OU model (see Note)}
  \item{r}{maximum-likelihood estimate of the slope parameter of early burst model}
  \item{z0}{maximum-likelihood estimate of \code{z0}, the value at the root of the tree}
  \item{nuisance}{maximum-likelihood estimate of \code{nuisance}, the unknown, nuisance contribution to measurement error (see details)}
  \item{convergence}{convergence diagnostics from \code{\link{optim}} function (see optim documentation)}
}

\details{

Note: if including known measurement error, the model fit incorporates this known error and, in addition, estimates an unknown, nuisance contribution to measurement error. The current implementation does not differentiate between the two, so, for instance, it is not possible to estimate the nuisance measurement error without providing the known, intraspecific error values.
 
}


\author{Jonathan Drury \href{mailto:jonathan.p.drury@gmail.com}{jonathan.p.drury@gmail.com}

		Julien Clavel 
}

\seealso{
\code{\link{fit_t_comp}}
\code{\link{sim_t_tworegime}}
}

\examples{
data(Cetacea)


\donttest{

regime<-c(rep("Odontocete",72),rep("Mystecite",15))
names(regime)<-c("Physeter_catodon","Kogia_breviceps","Kogia_simus","Platanista_gangetica","Platanista_minor","Tasmacetus_shepherdi","Berardius_arnuxii","Berardius_bairdii","Ziphius_cavirostris","Indopacetus_pacificus","Hyperoodon_ampullatus","Hyperoodon_planifrons","Mesoplodon_bidens","Mesoplodon_traversi","Mesoplodon_ginkgodens","Mesoplodon_europaeus","Mesoplodon_mirus","Mesoplodon_bowdoini","Mesoplodon_carlhubbsi","Mesoplodon_layardii","Mesoplodon_hectori","Mesoplodon_densirostris","Mesoplodon_stejnegeri","Mesoplodon_grayi","Mesoplodon_perrini","Mesoplodon_peruvianus","Lipotes_vexillifer","Inia_geoffrensis","Pontoporia_blainvillei","Delphinapterus_leucas","Monodon_monoceros","Neophocaena_phocaenoides","Phocoena_phocoena","Phocoena_dalli","Phocoena_sinus","Phocoena_dioptrica","Phocoena_spinipinnis","Orcinus_orca","Orcaella_brevirostris","Grampus_griseus","Pseudorca_crassidens","Feresa_attenuata","Peponocephala_electra","Globicephala_macrorhynchus","Globicephala_melas","Lagenorhynchus_albirostris","Lagenorhynchus_acutus","Lissodelphis_borealis","Lissodelphis_peronii","Cephalorhynchus_hectori","Cephalorhynchus_commersonii","Cephalorhynchus_eutropia","Lagenorhynchus_obscurus","Lagenorhynchus_obliquidens","Cephalorhynchus_heavisidii","Lagenorhynchus_australis","Lagenorhynchus_cruciger","Steno_bredanensis","Sotalia_fluviatilis","Sotalia_guianensis","Lagenodelphis_hosei","Stenella_longirostris","Stenella_attenuata","Tursiops_aduncus","Tursiops_truncatus","Sousa_chinensis","Stenella_clymene","Stenella_coeruleoalba","Stenella_frontalis","Delphinus_tropicalis","Delphinus_capensis","Delphinus_delphis","Balaena_mysticetus","Eubalaena_australis","Eubalaena_glacialis","Eubalaena_japonica","Caperea_marginata","Eschrichtius_robustus","Balaenoptera_acutorostrata","Balaenoptera_bonaerensis","Balaenoptera_physalus","Megaptera_novaeangliae","Balaenoptera_musculus","Balaenoptera_omurai","Balaenoptera_borealis","Balaenoptera_brydei","Balaenoptera_edeni")
regime.map<-phytools::make.simmap(Cetacea,regime)

data<-sim_t_tworegime(Cetacea,pars=c(sig2=0.01,r1=-0.01,r2=-0.02),root.value=0,Nsegments=1000,model="EB",regime.map=regime.map)
error<-rep(0.05,length(Cetacea$tip.label))
names(error)<-names(regime)

#Fit single-regime models
BM1.fit<-fit_t_standard(regime.map,data,model="BM",error,two.regime=FALSE)
OU1.fit<-fit_t_standard(regime.map,data,model="OU",error,two.regime=FALSE)
EB1.fit<-fit_t_standard(regime.map,data,model="EB",error,two.regime=FALSE)

#Now fit models that incorporate biogeography, NOTE these models take longer to fit
BM2.fit<-fit_t_standard(regime.map,data,model="BM",error,two.regime=TRUE)
OU2.fit<-fit_t_standard(regime.map,data,model="OU",error,two.regime=TRUE)
EB2.fit<-fit_t_standard(regime.map,data,model="EB",error,two.regime=TRUE)
}

}
