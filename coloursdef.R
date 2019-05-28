

WT=rgb(130/255,130/255,130/255)

APH=rgb(240/255,59/255,32/255)

CDC6=rgb(8/255,69/255,148/255)

cols3=c(WT, APH, CDC6)
names(cols3)=c('WT', 'APH', 'CDC6')


#CONST=rgb(35/255,139/255,69/255)
COMM='darkgreen'

APHR=rgb(254/255,178/255,76/255)

CDC6R=rgb(107/255,174/255,214/255)

APHCDC6R='green'

colsresp=c(COMM, APHR, CDC6R, APHCDC6R)
names(colsresp)=c('COMM', 'APH-R', 'CDC6-R', 'APH+CDC6-R')

colWTef='black'


colstot=c(cols3,colsresp)

sel=c('WT', 'APH', 'CDC6', 'APH-R', 'CDC6-R','APH+CDC6-R' , 'COMM',  'WT-nonCOMM','ALL-ORI')

cols=c(colstot,c( 'grey',  'black'))
names(cols)[8:9]=c( 'WT-nonCOMM', 'ALL-ORI')
cols
