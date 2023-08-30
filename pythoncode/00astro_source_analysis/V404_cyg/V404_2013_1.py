#importvla(archivefiles='AH823_A030729.xp1',vis='V404_2013_1.ms')

vis='V404_2013_1.ms'

flagdata(vis='V404_2013_1.ms',antenna='VA14',timerange='2003/07/29/00:16:15.0~14:13:57.0')
flagdata(vis='V404_2013_1.ms',antenna='VA03',timerange='2003/07/29/00:36:00.0~00:46:00.0')
flagdata(vis='V404_2013_1.ms',antenna='VA03',timerange='2003/07/29/01:03:40.0~01:08:00.0')
flagdata(vis='V404_2013_1.ms',antenna='VA03',timerange='2003/07/29/01:08:00.0~14:13:57.0')
flagdata(vis='V404_2013_1.ms',antenna='VA22',timerange='2003/07/29/12:59:30.0~13:17:33.0')
flagdata(vis='V404_2013_1.ms',antenna='VA16',timerange='2003/07/29/13:25:35.0~14:13:57.0')
flagdata(vis='V404_2013_1.ms',antenna='VA17',timerange='2003/07/29/14:01:00.0~14:13:57.0')
flagdata(vis='V404_2013_1.ms',mode='quack',quackinterval=10.0,quackmode='beg')
clearstat()


gencal(vis='V404_2013_1.ms',caltable='V404_2013_1.antpos',caltype='antpos')

setjy(vis='V404_2013_1.ms',field='1331+305',model='3C286_C.im',usescratch=False,scalebychan=True,spw='')







