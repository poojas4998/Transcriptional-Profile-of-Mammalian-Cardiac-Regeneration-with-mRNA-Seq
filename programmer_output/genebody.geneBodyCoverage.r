accepted_hits <- c(0.0,0.09524563085852311,0.1616940569598428,0.2035313716876168,0.24648980677994392,0.2752114303889752,0.30128535190110206,0.3305391927618454,0.35259078004468736,0.3772648757083764,0.4129809047378358,0.42938673771294594,0.4514091947615966,0.4618315201730808,0.49035631787541395,0.5148808258500124,0.5267927297233888,0.5355003078357181,0.5591194166395047,0.5805403737172856,0.6003961711850022,0.6167106769394044,0.6385268604378352,0.6494743173683903,0.6680436607099274,0.682957553312265,0.6937318034449258,0.7198419409238851,0.7233249721688167,0.7229659617149993,0.7130286153376194,0.7262757861620366,0.7344802772568452,0.7527071435207273,0.7685327337228912,0.7762176044240165,0.7816192946074425,0.8121446308254564,0.807208237085465,0.8136987681847455,0.8010633322783465,0.808580507306965,0.8200460100239497,0.8295456155848328,0.8459333405765267,0.8472079851483044,0.8502517009694857,0.8751383686124089,0.8916276557719591,0.9029569549615087,0.8997683752729976,0.8986724486245018,0.9004360087485179,0.9157600865404146,0.9121652581804783,0.9122605219193203,0.905670003259437,0.9076303892901513,0.9032403842671541,0.9199666498183691,0.9222049540293412,0.9279585689339437,0.9420394942991345,0.9415678994266855,0.9647358832523197,0.9793773687997871,0.9839374314061715,0.9834729222663637,0.969399082633814,0.9567392403149845,0.9593523010523101,0.9492472590024233,0.9630652249562652,0.9613906301420768,0.9519217294226388,0.9602238461671697,0.9691676153675369,0.9729576950269179,0.9886903834011201,1.0,0.994579414529531,0.9817455781091802,0.9728191295886024,0.970479263209774,0.9648098897932383,0.96498624580564,0.9423252855156603,0.9373684218813723,0.9335153579318478,0.9115834408002784,0.8936793690233813,0.895338217765034,0.8793024174945164,0.862042360084336,0.8251957630468021,0.7879736221792879,0.7460772596795044,0.6998215969981687,0.6144322202434658,0.43546314710453343)


pdf("/projectnb/bf528/users/group_5/project_2/programmer/output/genebody.geneBodyCoverage.curves.pdf")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(1)
plot(x,accepted_hits,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
dev.off()
