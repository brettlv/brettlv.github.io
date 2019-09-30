import os,sys


def get_file_list(pwdpath):
    dirname=os.listdir(pwdpath)
    dirname_select=[]
    for i in dirname:
        if i.endswith('.png') or i.endswith('.html'):
            dirname_select.append(i)
    return dirname_select




def main(argv):
    print(sys.argv)
    pwdpath=os.getcwd()
    print(pwdpath)
    print(pwdpath.split('/'))
    dirname_select=get_file_list(pwdpath)
    with open('Readme.md','w+') as f:
        f.write(pwdpath.split('/')[-2]+':'+pwdpath.split('/')[-1]+'\n')
        f.write('- [return back](../) \n filelist: \n')
        for i in dirname_select:
            f.write('- [%s](./%s)'%(i,i)+'\n')

if __name__ =="__main__":
    main(sys.argv)








