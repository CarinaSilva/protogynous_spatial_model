# protogynous_spatial_model

This code is used to implement a generic protogynous fish population model with spatial management (marine reserves + fishing)

Citation:
Easter EE, White JW (in review) Spatial management for protogynous sex-changing fishes: a general framework for coastal systems. Marine Ecology Progress Series

Language:
Matlab .m-files. Created and run with version 2015a

Contents:
# Lower-level functions called by other functions:
LifeHistory_Params.m: biological parameters for California sheephead
Spatial_Params.m: seascape parameters (reserve configuration, spatial array size, etc.)
Gonochore_F_FLEP.m: Calculates equivalent values of FLEP (fraction of unfished lifetime egg production) for a range of values of F, the fishing mortality rate
Find_F.m: Given a value of FLEP, returns the equivalent value of F (for a gonochore)
Spatial_Model.m: Implements the population dynamics model, given spatial and life history parameters

# Higher-level functions used to run the simulations in the paper:
Spatial_Struct.m: Runs the model for particular combinations of Phi and FLEP. Used to generate simulations in Fig. 2-3
Network_MinCR_Persist.m: Calculates the minimum fraction of reserves needed for population persistence (Fig. 4-5)
Self_MinRW_Persist.m: Calculates the minimum reserve width needed for population persistence (Fig. 4-5)

# Figure-generating code
PHI_Fert_fig.m: Makes Fig. 1
Network_SexRatio_Abundance_plots.m: Makes Fig. 2-3
Self_SexRatio_Abundance_plots.m: Makes equivalent versions of Fig. 2-3 for self-persistent case; not included in manuscript
Self_MinRW_Persist_fig.m: Makes self-persistence panels in Fig. 4-5
Network_MinCR_Persist_fig.m: Makes network-persistence panels in Fig. 4-5


