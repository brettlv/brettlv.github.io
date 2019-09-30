import os,sys



def get_file_list(pwdpath):
    dirname=os.listdir(pwdpath)
    dirname_select=[]
    dirs=[]
    for i in dirname:
        if i.endswith('.png') or i.endswith('.html'):
            dirname_select.append(i)
        if os.path.isdir(os.path.join(pwdpath,i)):
            dirs.append(i)

    return dirname_select,dirs



def main(argv):
    print(sys.argv)
    pwdpath=os.getcwd()
    print(pwdpath)
    print(pwdpath.split('/'))
    dirname_select,dirs=get_file_list(pwdpath)
    with open('Readme.md','w+') as f:
        f.write(pwdpath.split('/')[-2]+':'+pwdpath.split('/')[-1]+'\n')
        f.write('- [return back](../) \n dirslist: \n')
        
        for i in dirs:
            f.write('- [%s](./%s)'%(i,i)+'\n')
        
        f.write('\n fileslist: \n')
        for i in dirname_select:
            f.write('- [%s](./%s)'%(i,i)+'\n')

if __name__ =="__main__":
    main(sys.argv)








