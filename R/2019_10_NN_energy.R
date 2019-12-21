.encode <- function(reads)
{
  # if (nchar(reads)[1]!=seqLen) {reads=paste0(reads, paste0(rep("N",seqLen-nchar(reads)[1]),collapse=""))}
  seqLen=nchar(reads)[1]
  pos= reads %>% fjComm::seqFregments(1L)
  posA= ifelse(pos=="A",1,0)
  posC= ifelse(pos=="C",1,0)
  posG= ifelse(pos=="G",1,0)
  posT= ifelse(pos=="T",1,0)
  arr=array(0, dim=c(length(reads),seqLen,4))
  arr[,,1]=posA
  arr[,,2]=posC
  arr[,,3]=posG
  arr[,,4]=posT
  arr
}

.seed_encode<- function(seed, kernel_len=NULL)
{
  pos= seed %>% fjComm::seqFregments(1L)
  posA= ifelse(pos=="A",1,ifelse(pos=="N",0.25,0))
  posC= ifelse(pos=="C",1,ifelse(pos=="N",0.25,0))
  posG= ifelse(pos=="G",1,ifelse(pos=="N",0.25,0))
  posT= ifelse(pos=="T",1,ifelse(pos=="N",0.25,0))
  arr=array(0, dim=c(nchar(seed[1]),4))
  arr[,1]=posA
  arr[,2]=posC
  arr[,3]=posG
  arr[,4]=posT
  if(is.null(kernel_len)) {return(arr) }else
  {
    seedlen=nchar(seed)
    front_Ns=((kernel_len-seedlen)/2) %>% base::floor()
    end_Ns=((kernel_len-seedlen)/2) %>% base::ceiling()

    frontMat=matrix(0.25,nrow = front_Ns,ncol = 4)
    endMat=matrix(0.25,nrow = end_Ns,ncol = 4)
    return(rbind(frontMat,arr,endMat))
  }
}


NN_PWM_model_gen <- function(bseqs, ubseqs,  seqLen=30, both_strand=T,seed="AAAAAAA",TFname="Test", ker_size=14L)
{
  #seq_size=120000
  p_load(keras)

  if(both_strand){seqLen=seqLen*2}


  nrow_b= length(bseqs) %>% as.numeric()
  nrow_ub= length(ubseqs) %>% as.numeric()
  if(both_strand)
    {train_data= c(paste0(bseqs,bseqs %>% revComp()),paste0(ubseqs,ubseqs %>% revComp())) %>% .encode()}else
    {train_data= c(bseqs,ubseqs) %>% .encode()}

  train_label=c(rep(1,nrow_b),rep(0,nrow_ub))#cbind(c(rep(1,seq_size),rep(0,seq_size)), c(rep(0,seq_size),rep(1,seq_size)))

    set.seed(100)
    new_order=gtools::permute(1:length(train_label))
    train_label=train_label[new_order]
    train_data=train_data[new_order,,]


  bkend=backend()
  custom_exp= function(x){bkend$exp(x)}
  custom_log= function(x){bkend$log(x)}
  model2=keras_model_sequential() %>%
    layer_conv_1d(filter=1,kernel_size=ker_size, padding="valid",input_shape=list(seqLen,4),use_bias = F) %>%
    layer_activation(custom_exp) %>%
    layer_average_pooling_1d(pool_size = seqLen-ker_size+1) %>%
    layer_activation(custom_log) %>%
    layer_flatten() %>%
    layer_activation(activation = "sigmoid")


  weights=model2$get_weights();
  weights[[1]][,,1]=.seed_encode(seed,kernel_len =ker_size);
  model2$set_weights(weights)

  return(list(train_data=train_data, train_label=train_label,NNmodel=model2))
}






NN_PWM_default_train <- function(train_data, train_label, NNmodel, crude_train=TRUE, crude_epochs = 10, lr_crude=1, fine_train=TRUE, fine_epochs = 10, lr_fine=0.01, batch_size = 12000,validation_split=0.2)
{
  p_load(keras)

  if (crude_train)
  {
    #compiling the defined model with metric = accuracy and optimiser as adam.
    NNmodel %>% compile(
      loss = 'binary_crossentropy',
      optimizer = keras::optimizer_adam(lr=lr_crude),
      metrics = c(acc="accuracy")
    )

    #fitting the model on the training dataset
    history=NNmodel %>% fit(train_data, train_label, epochs = crude_epochs, batch_size = batch_size, validation_split=validation_split)
  }


  if (fine_train)
  {
    NNmodel %>% compile(
      loss = 'binary_crossentropy',
      optimizer = keras::optimizer_adam(lr=lr_fine),
      metrics = c(acc="accuracy")
    )

    #fitting the model on the training dataset
    history=NNmodel %>% fit(train_data, train_label, epochs = fine_epochs, batch_size = batch_size, validation_split=validation_split)
  }
  allweights2=get_weights(NNmodel);#allweights2[[1]][,,1] %>% t %>% sweep(2,colMins(.),FUN = "-") %>% {.+0.1} %>% fjComm::plotMotif_pfmMat()

  energy1=allweights2[[1]][,,1] %>% t %>% set_rownames(c("A","C","G","T")) %>% sweep(2,colMeans(.),FUN = "-") %>% `-`; prob1=exp(-energy1)
  # fjComm::plotMotif_pfmMat(prob1,name = "energy model",ic.scale = F)
  fjComm::plotMotif_pfmMat(prob1,name = "sigmoid(energy)--binary_crossEntropy ",ic.scale = T)

  return(prob1)
}



