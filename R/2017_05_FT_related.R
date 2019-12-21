# Welch window applied, minus mean to reduce signal around 0, normalized with average and length of the vector
# fit spline
# take 0.08-0.12 as signal
# take the range 0.14-0.3 /bp-1 as background to normalize
# check 10.2-bp periodicity
# default zero-padding factor of 0
# phase = avg(phase from 1/10.2 +- 0.0002)


win_Welch_coeff <- function(point_num)
{
  normFactor= 1-( (0:(point_num-1)-(point_num-1)/2) / ((point_num-1)/2) )**2;
  return (normFactor);
}

fft_raw <- function(Vect,win_Welch=TRUE,padding_ratio= 0)
{
  # average value (used to normalize power, as MI vibration scale with MI strength)
  Vect_avg= mean(Vect)
  Vect= Vect - Vect_avg

  if(win_Welch) Vect= win_Welch_coeff(length(Vect))*Vect
  Vect= c(Vect,rep(0,length(Vect)*padding_ratio))

  ft_result= fft(Vect)

  norm_pow= Mod(ft_result) / (Vect_avg * length(Vect)) #normalize with avg value, and with length of the ligand
  phase_angle= (ft_result %>% Arg())/pi  #in c(-1,1)

  # gen axis
  data_point_num=length(Vect); # points of data for FFT
  step=1; # step distance between data points

    sampling_rate= 1 #hz, so the upper limit of freq we can inspect is 0.5hz (upper bound of x-axis), which is 2bp periodicity
    max_x_axis= sampling_rate/2
    points_x_axis= (length(Vect)/2) %>% base::ceiling()
  x_axis= seq(0,0.5,length.out = points_x_axis)
  norm_pow= norm_pow[1:points_x_axis]
  phase_angle= phase_angle[1:points_x_axis]
  return(data_frame(x_axis=x_axis,norm_pow=norm_pow,phase_angle=phase_angle))
}


# calc sn ratio and phase angle at closest to 10bp
# phase yet to be extrapolated between 9.5 and 10.6 !!!!!!!!!!!!!!!!!!!!!
# take 0.08-0.12 as signal
# take the range 0.14-0.3 /bp-1 as background to normalize
## corr by substract ref region##
fft_snRatio_phase <- function(Vect, periodicity=10.2, padding_ratio= 0, phase_int_width=0.0002)
{
  fft_result=fft_raw(Vect,win_Welch = T,padding_ratio = padding_ratio)
#   ggplot()+geom_line(data = fft_result,aes(x_axis,norm_pow)) + geom_line(data = fft_raw(Vect,win_Welch = F),aes(x_axis,norm_pow),color="red")

  pacman::p_load(MESS)
  sig= MESS::auc(x = fft_result$x_axis, y = fft_result$norm_pow, from = 0.08, to = 0.12,type = "linear")
  bk_density= MESS::auc(x = fft_result$x_axis, y = fft_result$norm_pow, from = 0.14, to = 0.3,type = "linear") / (0.3-0.14)
    # bk_density_near=  MESS::auc(x = fft_result$x_axis, y = fft_result$norm_pow, from = 0.04, to = 0.08,type = "linear") / (0.08-0.04)

  # sig_bk_ratio= sig/bk
  sig_sub_bk= sig-bk_density*(0.12-0.08)
    ind_phase= with(fft_result,(x_axis-1/periodicity) %>% abs() %>% which.min())
  phase_angle_in_pi= fft_result$phase_angle[ind_phase]
  freq_of_phase_angle= fft_result$x_axis[ind_phase]
  # phase can not take segment avg because abrupt change
  return(c(sig_bk_ratio=sig_sub_bk, phase_angle_in_pi=phase_angle_in_pi,freq_of_phase_angle=freq_of_phase_angle))
}

## corr by baseline sub
fft_snRatio_phase_bl_sub <- function(Vect, periodicity=10.2, padding_ratio= 0, phase_int_width=0.0002,baseline_hwm=10, baseline_start_freq=0.04)
{
  fft_result=fft_raw(Vect,win_Welch = T,padding_ratio = padding_ratio) %>% dplyr::filter(x_axis >=baseline_start_freq)
  # ggplot()+geom_line(data = fft_result,aes(x_axis,norm_pow)) + geom_line(data = fft_raw(Vect,win_Welch = F),aes(x_axis,norm_pow),color="red")

  pacman::p_load(baseline)
  bl= baseline.medianWindow(matrix(fft_result$norm_pow,nrow = 1),hwm = baseline_hwm)$baseline[1,] #$corrected[1,]
  pacman::p_load(MESS)
  sig= MESS::auc(x = fft_result$x_axis, y = fft_result$norm_pow -bl, from = 0.08, to = 0.12,type = "linear")
  # bk= MESS::auc(x = fft_result$x_axis, y = fft_result$norm_pow, from = 0.14, to = 0.3,type = "linear")


  # sig_bk_ratio= sig/bk
  sig_bk_ratio= sig
  ind_phase= with(fft_result,(x_axis-1/periodicity) %>% abs() %>% which.min())
  phase_angle_in_pi= fft_result$phase_angle[ind_phase]
  freq_of_phase_angle= fft_result$x_axis[ind_phase]
     # phase can not take segment avg because abrupt change
  # browser()
  # return(ggplot(fft_result,aes(x_axis,norm_pow))+geom_line()+geom_line(aes(y=norm_pow-bl),color="red"))
  return(c(sig_bk_ratio=sig_bk_ratio, phase_angle_in_pi=phase_angle_in_pi,freq_of_phase_angle=freq_of_phase_angle))
}






# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/P24_Trulig147v1IIIFJ4-5-CAP2r-147c4-P24-SHOXIIIc4-GGTATGAA_S384_L003_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/K20_Trulig147v1IIIFJ4-5-CAP2r-147c4-K20-PROX1IIIc4-ACTTGTAG_S260_L003_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"

# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/E4_Trulig147v1IIIFJ4-5-CAP2r-147c4-E4-PROX2IIIc4-CCGGCTTT_S100_L003_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"
# file="~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/F10_Trulig147v1IIIFJ4-5-CAP2r-147c4-F10-NOTFIIIc4-TAGATACT_S130_L003_R1_001.peared_trimmed.fq.gz_3mer_foldchn_maxBias.csv"

# tfmi=rio::import(file)
# fjComm::gg_heat2D_diag(tfmi)
# tfmi %<>% dplyr::filter(pos2-pos1==3) %>% .$topMIsum
# fft_snRatio_phase(tfmi) %>% print()
# fft_snRatio_phase(runif(100),baseline_hwm = 20) %>% print()
# print(fft_snRatio_phase(c(tfmi)))
