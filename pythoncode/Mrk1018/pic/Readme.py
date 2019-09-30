import os,sys



def get_file_list(pwdpath):
    dirname=os.listdir(pwdpath)
    dirname_png=[]
    dirname_html=[]
    dirname_ipynb=[]
    dirs=[]
    
    for i in dirname:
        if i.endswith('.png'):
            dirname_png.append(i)
        if i.endswith('.html'):
            dirname_html.append(i)
        if i.endswith('.ipynb'):
            dirname_ipynb.append(i)
        if os.path.isdir(os.path.join(pwdpath,i)):
            dirs.append(i)
    return dirs,dirname_png,dirname_html,dirname_ipynb



def main(argv):
    print(sys.argv)
    pwdpath=os.getcwd()
    print(pwdpath)
    print(pwdpath.split('/'))
    dirs,dirname_png,dirname_html,dirname_ipynb=get_file_list(pwdpath)
    with open('Readme.md','w+') as f:
        f.write(pwdpath.split('/')[-2]+':'+pwdpath.split('/')[-1]+'\n')
        f.write('- [return back](../) \n')
        f.write('\n dirslist: \n')
        for i in dirs:
            f.write('- [%s](./%s)'%(i,i)+'\n')
        
        f.write('\n files_list: \n')
        f.write('\n ipynb_list: \n')
        for i in dirname_ipynb:
            f.write('- [%s](./%s)'%(i,i)+'\n')
        f.write('\n html_list: \n')
        for i in dirname_html:
            f.write('- [%s](./%s)'%(i,i)+'\n')
        f.write('\n png_list: \n')
        for i in dirname_png:
            f.write('- [%s](./%s)'%(i,i)+'\n')

if __name__ =="__main__":
    main(sys.argv)








