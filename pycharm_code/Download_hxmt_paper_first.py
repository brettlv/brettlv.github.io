# coding = UTF-8
# 爬取html链接中对应的PDF文档,
# 网址：https://blogs.cornell.edu/arxiv/2017/10/16/gw170817/
# hxmt: https://www.sciencedirect.com/journal/journal-of-high-energy-astrophysics/special-issue/10B0FLKPH6C

from urllib.request import urlretrieve
from urllib.request import urlopen
from bs4 import BeautifulSoup 
import numpy as np
import pandas as pd
import urllib
import re
import os
from lxml.html import parse
import urllib.request
import re
import os
import requests,time
# open the url and read


def getHtml(url):
    page = urllib.request.urlopen(url)
    html = page.read()
    page.close()
    return html
# compile the regular expressions and find
# all stuff we need
def getUrl(html):
    reg = r'(?:href|HREF)="?((?:https://)?.+?\.pdf)' #匹配了
    url_re = re.compile(reg)
    url_lst = url_re.findall(html.decode('UTF-8')) #返回匹配的数组
    return(url_lst)


#raw_url = 'https://science.nrao.edu/science/meetings/2018/16th-synthesis-imaging-workshop/talks/'
#print(url_lst)
#path='/Users/brettlv/16th_synthesis_imaging_workshop/'


url_hxmt_first='https://www.sciencedirect.com/journal/journal-of-high-energy-astrophysics/special-issue/10B0FLKPH6C'
#raw_url = 'https://science.nrao.edu/science/meetings/2016/vla-data-reduction/program/'
raw_url=url_hxmt_first
path='/Users/brettlv/Downloads/2020_1th_HXMT_papers/'

isExists = os.path.exists(path)
if not isExists:
    os.makedirs(path)
os.chdir(path)
print(os.getcwd())

html = getHtml(raw_url)
url_lst = getUrl(html)
#print(url_lst)

def downpdf(pdfUrl):
    filename = pdfUrl.split('/')[-1]
    try:
        res = requests.get(pdfUrl,stream=True,timeout=10)
        if str(res.status_code)[0] == "4":
            print(str(res.status_code), ":" , pdfUrl)
            return False
    except Exception as e:
        print("抛出异常：", pdfUrl)
        print(e)
        return False
    #block_sz = 8192*16
    with open(filename, "wb") as f:
        for chunk in res.iter_content(chunk_size=102400):
            if chunk:
                f.write(chunk)
        print(filename,'downloaded')        
        #f.write(res.content)    
    return True

for url in url_lst[:]:
    print(url)
    downpdf(url)