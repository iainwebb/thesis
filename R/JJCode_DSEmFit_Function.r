#####
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
#####
##  -> Method derived/recommended to JJ by Dr David Sexton at UK Met Office (Hence "DSEm")
##
##  -> This code produces an emulator fit with a reduced mean function form that includes only 
##      key parameters (removes some over-fitting) and also initiates the parameter estimation 
##      of the covariance function parameters in the emulator fit using values closer to the  
##      best MLE estimates, and therefore get better convergence to obtain a more robust and 
##      stable emulator model fit...
##
## Function Inputs:
##  -> "TrInputs" is a table containing the training input combinations to fit the emulator 
##      with (1 column per input, number of rows = number input combs).
##  -> "TrOutput" is a vector containing the corresponing model output for the training runs 
##       (length = number of training runs; values in order of training runs) 
##  -> "InputNames" is a vector of the input names. eg: InputNames=c("Name1",...,"NameN")...  
##       -> This should be the same as the output from: colnames(TrInputs) 
##       -> The table TrInputs should have column headers the same as "InputNames", and so if 
##           this isn't the case, need to do: colnames(TrInputs)<-InputNames prior to running.
##
## To RUN THIS CODE:
##  -> Source this code file in the R workspace using the command:
##       source("FILEPATH_TO_WHERE_YOU_HAVE_SAVED_THIS FILE/JJCode_DSEmFit_Function.r")
##  -> Then, run the code as a function call: 
## EmulatorModel<-JJCode_FitDSEm_EmModel(TrInputs=YOUR_INPUT_TABLE,TrOutput=YOUR_OUTPUT_VECTOR,InputNames=YOUR_INPUT_NAMES_VECTOR)
#####
###############################################################################################
#####
JJCode_FitDSEm_EmModel<-function(TrInputs,TrOutput,InputNames,nuggetIn=NULL,nuggetEstimIn=FALSE,noiseVarIn=NULL,DKcontrol_list=list(trace=3, maxit=500, REPORT=1, factr=1e7, pgtol=0.0, pop.size=100))
{
#####
## -> Fit initial emulator with const mean & save covariance parameter estimates...
#####
    EmModInit<-km(~1,design=data.frame(TrInputs),response=TrOutput,covtype='matern5_2',optim.method='BFGS',nugget=nuggetIn,nugget.estim=nuggetEstimIn,noise.var=noiseVarIn,control=DKcontrol_list)
    coefs0<-EmModInit@covariance@range.val
#####
## -> Make a form of the input/output data for the lm() regression function...
## -> Make the linear regression formula for fitting, regressing on all inputs...
## -> Fit this linear regression model:
#####
    RegrData<-data.frame(y=TrOutput,x=TrInputs)
    colnames(RegrData)<-c("y",InputNames)
    RegrStartForm<-as.formula(paste("y ~ ", paste(InputNames, collapse= "+")))
    Startlm<-lm(RegrStartForm,data=RegrData)
#####
## -> Reduce this linear regression model down using step() with BIC criterion to only regress
##     the output on the key input parameters it is sensitive to...
## -> Make the final formula for fitting the final Em Model with...
#####
    nTr<-dim(TrInputs)[1]
    Steplm<-step(Startlm, direction="both", k=log(nTr), trace=FALSE)
    Labels<-labels(terms(Steplm))
    Labels<-Labels[!(Labels %in% c('y','response'))]
    if(length(Labels) > 0){
      RegrEndform<-as.formula(paste("~ ", paste(Labels, collapse= "+")))
    }else{
      RegrEndform<-as.formula("~ 1")
    }
#####
## Fit the final form of the emulator model, using the reduced set of key inputs for the mean 
##   function & initialising the covariance pars with the initial cov par values, "coefs0"...
#####
    EmMod<-km(RegrEndform,design=data.frame(TrInputs),response=TrOutput,covtype='matern5_2',optim.method='BFGS',nugget=nuggetIn,nugget.estim=nuggetEstimIn,noise.var=noiseVarIn,parinit=coefs0,control=DKcontrol_list)
#####
## Return the emulator:
#####
  return(EmMod)
}
#####
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
#####

