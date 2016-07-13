##Run t.test on deviations from expectations

#Time subsets
T1_diff=subset(data,data$Time==1)
T2_diff=subset(data,data$Time==2)
T1_A_Dark_diff=subset(diff,diff$Clade=="A"|diff$Env=="Dark"|diff$Time==1)
T1_B_Dark_diff=subset(diff,diff$Clade=="B"|diff$Env=="Dark"|diff$Time==1)
T1_D_Dark_diff=subset(diff,diff$Clade=="D"|diff$Env=="Dark"|diff$Time==1)
T2_A_Dark_diff=subset(diff,diff$Clade=="A"|diff$Env=="Dark"|diff$Time==2)
T2_B_Dark_diff=subset(diff,diff$Clade=="B"|diff$Env=="Dark"|diff$Time==2)
T2_D_Dark_diff=subset(diff,diff$Clade=="D"|diff$Env=="Dark"|diff$Time==2)
T1_A_Light_diff=subset(diff,diff$Clade=="A"|diff$Env=="Light"|diff$Time==1)
T1_B_Light_diff=subset(diff,diff$Clade=="B"|diff$Env=="Light"|diff$Time==1)
T1_D_Light_diff=subset(diff,diff$Clade=="D"|diff$Env=="Light"|diff$Time==1)
T2_A_Light_diff=subset(diff,diff$Clade=="A"|diff$Env=="Light"|diff$Time==2)
T2_B_Light_diff=subset(diff,diff$Clade=="B"|diff$Env=="Light"|diff$Time==2)
T2_D_Light_diff=subset(diff,diff$Clade=="D"|diff$Env=="Light"|diff$Time==2)