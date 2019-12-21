# arrange a vector to 384 plate format
excel_vect_to_384 <- function(arranged_384_vect)
{
  arranged_384_vect %>% matrix(ncol = 24,byrow=T) %>% {cbind(LETTERS[1:16],.)} %>% {rbind(0:24,.)}
}

excel_384region_to_list<- function() readClipBoard() %>% as.matrix() %>% t %>% as.character() %>% {data.frame(well=fjComm::gen_384_label(),TF=.)}
excel_96region_to_list<- function() readClipBoard() %>% as.matrix() %>% t %>% as.character() %>% {data.frame(well=fjComm::gen_96_label(),TF=.)}


# when combining 4x96 to form 384, convert 96 well index to 384 well index
wellInd_96_to_384=Vectorize(FUN =
                              function (Ind96="A1", quad=4)
                              {
                                row96= stringr::str_sub(Ind96, 1,1)
                                row96 = utf8ToInt(row96) - 64 #ASC("A")=65
                                col96 = stringr::str_sub(Ind96, 2, -1) %>% as.integer()
                                if(quad==1)
                                {
                                  row384 = (row96 - 1) * 2 + 1
                                  col384 = (col96 - 1) * 2 + 1
                                }
                                else if(quad==2)
                                {
                                  row384 = (row96 - 1) * 2 + 1
                                  col384 = (col96) * 2
                                }
                                else if(quad==3)
                                {
                                  row384 = (row96) * 2
                                  col384 = (col96 - 1) * 2 + 1
                                }
                                else if(quad==4)
                                {
                                  row384 = (row96) * 2
                                  col384 = (col96) * 2
                                }else{stop("quad should be 1-4")}

                                row384 = intToUtf8(64 + row384) #'change to letter
                                wellInd96To384 = paste0(row384,col384)
                                return(wellInd96To384)
                              }

                              ,vectorize.args = c("Ind96","quad") )




# arragne 4x96 plate wells into 1 384 plate well
# fjComm::clear_()
# letter2num <- function(x) { sapply(x, function(x) {utf8ToInt(x) - utf8ToInt("A") + 1L}) }
# num2letter <- function(x) {  sapply(x, function(x) {intToUtf8 (x + utf8ToInt("A") - 1L)})  }
#
# pacman::p_load(readxl)
# mainTable=read_excel("Information_for_pETG20A_PSC_FLTF1-4.xlsx", sheet = 1,skip = 1) %>% mutate(well384="") %>%
#   mutate(Plate=Plate %>% as.integer(),rows=Well %>% str_sub(1,1) %>% letter2num, cols=Well %>% str_sub(2,3) %>% as.integer()) %>%
#   mutate(row384=case_when(
#     Plate==1 | Plate==2 ~ 1+2*(rows-1),
#     Plate==3 | Plate==4 ~ 2*rows,
#     TRUE ~ 0# stop("not correct rows")
#   )) %>%
#   mutate(row384=num2letter(row384)) %>%
#   mutate(col384=case_when(
#     Plate==1 | Plate==3 ~ 1+2*(cols-1),
#     Plate==2 | Plate==4 ~ 2*cols,
#     TRUE ~ 0# stop("not correct rows")
#   ) %>% as.character()) %>%
#   mutate(col384=ifelse(nchar(col384)==1,paste0("0",col384),col384)) %>%
#   mutate(well384=paste0(row384,col384))
#
# fjComm::writeClipboard(mainTable)
