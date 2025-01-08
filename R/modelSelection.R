setGeneric(
    name="modelSelection",
    def=function(object="PhenotypicModel", data="numeric"){standardGeneric("modelSelection")}
)

setMethod(
    f="modelSelection",
    signature="PhenotypicModel",
    definition=function(object, data){
        message("*** Model selection with tip trait data ***\n")
        message("For each model in \"object\", fits the model and returns its AIC value in a recap table...\n")
        message("**WARNING** : This function relies on the standard R optimizer \"optim\".\nIt may not always converge well.\nPlease double check the convergence by trying\n distinct parameter sets for the initialisation.")

        aic <- c()
        names <- c()
        for(model in object){
            fit <- fitTipData(model, data)
            aic <- c(aic, 2*length(model@params0)+2*fit$value )
            names <- c(names, model@name)
        }
        names(aic) <- names

        return(sort(aic))
    }
)
