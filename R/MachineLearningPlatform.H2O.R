#' Build ensemble model based on h2o platform
#'
#' @param train.df row is sample, column is feature
#' @param train.group a vector
#' @param val.df Default NULL, not to evaluate in validation dataset. row is sample, column is feature
#' @param val.group a vector
#' @param candidates if not specify candidates feature, all column names will be used
#' @param nfolds Default 5. Used in internal model construction
#' @param seed Default 1
#' @param type Description for this run. Default "run"
#' @param metalearner_algorithm Default glm. Could be "AUTO", "deeplearning", "drf", "gbm", "glm", "naivebayes", "xgboost". If AUTO,  (GLM with non negative weights; if validation_frame is present, a lambda search is performed
#' @param glm TRUE. Default included glm
#' @param xg TRUE. Default included xgboost
#' @param dl TRUE. Default included deep learning
#' @param nb TRUE. Default included NaiveBayes
#' @param gbm TRUE. Default included gbm
#' @param rf TRUE. Default included random forests
#'
#' @return A list include fits, training risk score, aucs
#' @export
#'
#' @examples
buildEnsembleModel <- function(train.df = NULL, train.group = NULL,
                          val.df = NULL, val.group = NULL,
                          candidates = NULL, nfolds = 5,
                          seed = 1, type = "run", metalearner_algorithm = "glm",
                          glm = TRUE, xg = TRUE, dl = TRUE,
                          nb = TRUE, gbm = TRUE, rf = TRUE){

  if(!require("h2o")){
    BiocManager::install("h2o")
    require("h2o")
  }
  h2o.init(nthreads = -1, max_mem_size = "80G") #-1表示使用你机器上所有的核


  if(is.null(train.df) | is.null(train.group)){
    stop("Train.d and train.group is mandatory")
  }
  train.samples = rownames(train.df)
  train.df$Class = as.factor(train.group)
  train.df = as.h2o(train.df)

  # if has validation data set
  has.val = FALSE
  val.samples = NULL
  if(!is.null(val.df)&!is.null(val.group)){
    val.samples = rownames(val.df)

    val.df$Class = as.factor(val.group)
    val.df = as.h2o(val.df)
    has.val = TRUE
    warning("Have validation dataset.")
  }else{
    warning("No validation dataset.")
  }

  # if not specify features
  if(is.null(candidates)){
    candidates = colnames(train.df)
  }

  #
  base_models4ensemble = list()

  train.aucs = list()
  val.aucs = list()
  models = list()
  perf_list = list()

  cat("Number of training samples", length(train.samples),"\n")
  cat("Number of validation samples", length(val.samples),"\n")

  train.risk.scores = data.frame(Samples = train.samples,
                                 Class = train.group)
  if(has.val){
    val.risk.scores = data.frame(Samples = val.samples,
                                 Class = val.group)
  }


  ############# glm
  if(glm){

    glm_fit <- h2o.glm(x = candidates,
                       y = "Class",
                       training_frame = train.df,
                       model_id = paste0(type,".glm_fit"),
                       family = "binomial",
                       seed = seed,
                       keep_cross_validation_predictions = T,
                       lambda_search = F,
                       nfolds = nfolds)

    base_models4ensemble <- append(base_models4ensemble, glm_fit@model_id)

    models = append(models, list("glm"=glm_fit))

    glm_perf <- h2o.performance(model = glm_fit, newdata = train.df)
    perf_list = append(perf_list, list("train.glm" = glm_perf) )

    train.aucs <- append(train.aucs, list("glm" = h2o.auc(glm_perf)) )
    train.risk.scores$glm = unlist(as.data.frame(h2o.predict(glm_fit, train.df))[,c(3)])

    if(has.val){
      glm_perf <- h2o.performance(model = glm_fit, newdata = val.df)
      perf_list = append(perf_list, list("val.glm" = glm_perf) )

      val.aucs <- append(val.aucs, list("glm" = h2o.auc(glm_perf)) )
      val.risk.scores$glm = unlist(as.data.frame(h2o.predict(glm_fit, val.df))[,c(3)])
    }

  }



  ############# RF
  if(rf){

    rf_fit <- h2o.randomForest(x = candidates,
                               y = "Class",
                               training_frame = train.df,
                               model_id = paste0(type,".rf_fit"),
                               keep_cross_validation_predictions = T,
                               seed = seed,
                               nfolds = nfolds)

    base_models4ensemble <- append(base_models4ensemble, rf_fit@model_id)

    models = append(models, list("rf" = rf_fit))

    rf_perf <- h2o.performance(model = rf_fit, newdata = train.df)
    perf_list = append(perf_list, list("train.rf" = rf_perf) )

    train.aucs <- append(train.aucs, list("rf" = h2o.auc(rf_perf)) )
    train.risk.scores$rf = unlist(as.data.frame(h2o.predict(rf_fit, train.df))[,c(3)])

    if(has.val){
      rf_perf <- h2o.performance(model = rf_fit, newdata = val.df)
      perf_list = append(perf_list, list("val.rf" = rf_perf) )

      val.aucs <- append(val.aucs, list("rf" = h2o.auc(rf_perf)) )
      val.risk.scores$rf = unlist(as.data.frame(h2o.predict(rf_fit, val.df))[,c(3)])

    }

  }


  ########### GBM
  if(gbm){

    gbm_fit <- h2o.gbm(x = candidates,
                       y = "Class",
                       training_frame = train.df,
                       model_id = paste0(type,".gbm_fit"),
                       keep_cross_validation_predictions = T,
                       seed = seed,
                       nfolds = nfolds)

    base_models4ensemble <- append(base_models4ensemble, gbm_fit@model_id)

    models = append(models, list("gbm" = gbm_fit))

    gbm_perf <- h2o.performance(model = gbm_fit, newdata = train.df)
    perf_list = append(perf_list, list("train.gbm" = gbm_perf) )

    train.aucs <- append(train.aucs, list("gbm" = h2o.auc(gbm_perf)) )
    train.risk.scores$gbm = unlist(as.data.frame(h2o.predict(gbm_fit, train.df))[,c(3)])

    if(has.val){
      gbm_perf <- h2o.performance(model = gbm_fit, newdata = val.df)
      perf_list = append(perf_list, list("val.gbm" = gbm_perf) )

      val.aucs <- append(val.aucs, list("gbm" = h2o.auc(gbm_perf)) )
      val.risk.scores$gbm = unlist(as.data.frame(h2o.predict(gbm_fit, val.df))[,c(3)])

    }

  }


  ########### deep learning
  if(dl){

    dl_fit <- h2o.deeplearning(x = candidates,
                               y = "Class",
                               training_frame = train.df,
                               model_id = paste0(type,".dl_fit"),
                               epochs = 10,
                               keep_cross_validation_predictions = T,
                               nfolds = nfolds,
                               seed = seed)
    base_models4ensemble <- append(base_models4ensemble, dl_fit@model_id)

    models = append(models, list("dl" = dl_fit))

    dl_perf <- h2o.performance(model = dl_fit, newdata = train.df)
    perf_list = append(perf_list, list("train.dl" = dl_perf) )

    train.aucs <- append(train.aucs, list("dl" = h2o.auc(dl_perf)) )
    train.risk.scores$dl = unlist(as.data.frame(h2o.predict(dl_fit, train.df))[,c(3)])

    if(has.val){
      dl_perf <- h2o.performance(model = dl_fit, newdata = val.df)
      perf_list = append(perf_list, list("val.dl" = dl_perf) )

      val.aucs <- append(val.aucs, list("dl" = h2o.auc(dl_perf)) )
      val.risk.scores$dl = unlist(as.data.frame(h2o.predict(dl_fit, val.df))[,c(3)])

    }

  }

  ########## native bayes
  if(nb){

    nb_fit <- h2o.naiveBayes(x = candidates,
                             y = "Class",
                             training_frame = train.df,
                             model_id = paste0(type,".nb_fit"),
                             keep_cross_validation_predictions = T,
                             nfolds = nfolds,
                             seed = seed)
    base_models4ensemble <- append(base_models4ensemble, nb_fit@model_id)

    models = append(models, list("nb" = nb_fit))

    nb_perf <- h2o.performance(model = nb_fit, newdata = train.df)
    perf_list = append(perf_list, list("train.nb" = nb_perf) )

    train.aucs <- append(train.aucs, list("nb" = h2o.auc(nb_perf)) )
    train.risk.scores$nb = unlist(as.data.frame(h2o.predict(nb_fit, train.df))[,c(3)])

    if(has.val){
      nb_perf <- h2o.performance(model = nb_fit, newdata = val.df)
      perf_list = append(perf_list, list("val.nb" = nb_perf) )

      val.aucs <- append(val.aucs, list("nb" = h2o.auc(nb_perf)) )
      val.risk.scores$nb = unlist(as.data.frame(h2o.predict(nb_fit, val.df))[,c(3)])

    }

  }

  ############ XGBoots
  if(xg){

    xg_fit <- h2o.xgboost(x = candidates,
                          y = "Class",
                          training_frame = train.df,
                          model_id = paste0(type,".xg_fit"),
                          keep_cross_validation_predictions = T,
                          nfolds = nfolds,
                          seed = seed)
    base_models4ensemble <- append(base_models4ensemble, xg_fit@model_id)

    models = append(models, list("xg" = xg_fit))

    xg_perf <- h2o.performance(model = xg_fit, newdata = train.df)
    perf_list = append(perf_list, list("train.xg" = xg_perf) )

    train.aucs <- append(train.aucs, list("xg" = h2o.auc(xg_perf)) )
    train.risk.scores$xg = unlist(as.data.frame(h2o.predict(xg_fit, train.df))[,c(3)])

    if(has.val){
      xg_perf <- h2o.performance(model = xg_fit, newdata = val.df)
      perf_list = append(perf_list, list("val.xg" = xg_perf) )

      val.aucs <- append(val.aucs, list("xg" = h2o.auc(xg_perf)) )
      val.risk.scores$xg = unlist(as.data.frame(h2o.predict(xg_fit, val.df))[,c(3)])

    }


  }



  ############
  ensemble_fit <- h2o.stackedEnsemble(x = candidates,
                                      y = "Class",
                                      training_frame = train.df,
                                      model_id = paste0(type,".ensemble"),
                                      base_models = base_models4ensemble,
                                      seed = seed,
                                      metalearner_algorithm = metalearner_algorithm)

  models = append(models, list("ensemble" = ensemble_fit))

  ensemble_perf <- h2o.performance(model = ensemble_fit, newdata = train.df)
  perf_list = append(perf_list, list("train.ensemble" = ensemble_perf) )

  train.aucs <- append(train.aucs, list("ensemble" = h2o.auc(ensemble_perf)) )
  train.risk.scores$ensemble = unlist(as.data.frame(h2o.predict(ensemble_fit, train.df))[,c(3)])

  if(has.val){
    ensemble_perf <- h2o.performance(model = ensemble_fit, newdata = val.df)
    perf_list = append(perf_list, list("val.ensemble" = ensemble_perf) )

    val.aucs <- append(val.aucs, list("ensemble" = h2o.auc(ensemble_perf)) )
    val.risk.scores$ensemble = unlist(as.data.frame(h2o.predict(ensemble_fit, val.df))[,c(3)])

  }

  train.aucs = unlist(train.aucs)
  val.aucs = unlist(val.aucs)

  if(type != "run"){
    colnames(train.risk.scores)[3:ncol(train.risk.scores)] = paste0(type,".",colnames(train.risk.scores)[3:ncol(train.risk.scores)])
    if(has.val){
      colnames(val.risk.scores)[3:ncol(val.risk.scores)] = paste0(type,".",colnames(val.risk.scores)[3:ncol(val.risk.scores)])
    }
  }


  res = list(fits = models, seed = seed,
             train.aucs = train.aucs, val.aucs = val.aucs,
             performance.objects = perf_list,
             train.risk = train.risk.scores, val.risk = val.risk.scores)
  res

}
