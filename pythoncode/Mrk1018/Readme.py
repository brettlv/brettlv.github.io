import os,sys


def get_file_list(pwdpath):
    dirname=os.listdir(pwdpath)
    dirname_png=[]
    dirname_html=[]
    dirname_ipynb=[]
    dirname_pdf=[]
    dirname_pptx=[]
    dirs=[]
    for i in dirname:
        if i.endswith('.png'):
            dirname_png.append(i)
        if i.endswith('.html'):
            dirname_html.append(i)
        if i.endswith('.ipynb'):
            dirname_ipynb.append(i)
        if i.endswith('.pdf'):
            dirname_pdf.append(i)
        if i.endswith('.pptx'):
            dirname_pptx.append(i)
        if os.path.isdir(os.path.join(pwdpath,i)):
            dirs.append(i)
    return dirs,dirname_png,dirname_html,dirname_ipynb,dirname_pdf,dirname_pptx


'''<div align="center"><iframe src="https://view.officeapps.live.com/op/view.aspx?src=https://brettlv.github.io/radiolearnnote/v404_notebook_lvbing_20180104.pptx" frameborder="0" width="900" height="550" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe></div>'''

def main(argv):
    print(sys.argv)
    pwdpath=os.getcwd()
    print(pwdpath)
    print(pwdpath.split('/'))
    web=''
    for i in pwdpath.split('/')[5:]:
        web=web+i+'/'
    print(web)
    
    dirs,dirname_png,dirname_html,dirname_ipynb,dirname_pdf,dirname_pptx=get_file_list(pwdpath)
    if os.path.isfile(os.path.join(pwdpath,'Readme.md')):
        mdname='Readme_list.md'
    else:
        mdname='Readme.md'
    with open(mdname,'w+') as f:
        f.write(pwdpath.split('/')[-2]+':'+pwdpath.split('/')[-1]+'\n')
        f.write('- [return back](../) \n')
        f.write('\n dirslist: \n')
        for i in dirs:
            f.write('- [%s](./%s)'%(i,i)+'\n')
        
        f.write('\n files_list: \n')
        if len(dirname_ipynb)>0:
            f.write('\n ipynb_list: \n')
        for i in dirname_ipynb:
            f.write('- [%s](http://nbviewer.jupyter.org/github/brettlv/brettlv.github.io/tree/master/%s%s)'%(i,web,i)+'\n')
        if len(dirname_pdf)>0:
            f.write('\n pdf_list: \n')
        for i in dirname_pdf:
            f.write('- [%s](./%s)'%(i,i)+'\n')
        if len(dirname_pptx)>0:
            f.write('\n pptx_list: \n')
        for i in dirname_pptx:
            f.write('- [%s](https://view.officeapps.live.com/op/view.aspx?src=https://brettlv.github.io/%s%s)'%(i,web,i)+'\n')
        if len(dirname_html)>0:
            f.write('\n html_list: \n')
        for i in dirname_html:
            f.write('- [%s](./%s)'%(i,i)+'\n')
        if len(dirname_png)>0:
            f.write('\n png_list: \n')
        for i in dirname_png:
            f.write('- [%s](./%s)'%(i,i)+'\n')






if __name__ =="__main__":
    main(sys.argv)








