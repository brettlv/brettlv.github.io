import urllib.parse
import requests
from pyquery import PyQuery as pq
from urllib.parse import urlencode

headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.77 Safari/537.36',
    #'referer': 'https://so.toutiao.com/search?keyword=%E8%A1%97%E6%8B%8D&pd=atlas&dvpf=pc&aid=4916&page_num=0&search_json={%22from_search_id%22:%22202106100003510102121720341003A4ED%22,%22origin_keyword%22:%22%E8%A1%97%E6%8B%8D%22,%22image_keyword%22:%22%E8%A1%97%E6%8B%8D%22}',
    #'cookie': 'tt_webid=6970268611437856264; _S_DPR=1; _S_IPAD=0; MONITOR_WEB_ID=6970268611437856264; ttwid=1%7CobYZC0vIyDzWzgZETeBnLXiczCSfyj93HRpnKzDmjv0%7C1623254627%7C0320a1f1cb08676972a90f718fbd310044d6554255abcf4aa2e226f2262a7a9f; _S_WIN_WH=171_657',
}

#获取首页
def get_page(page_num):
    global headers
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.77 Safari/537.36',
    }
    params = {
        'keyword':urllib.parse.unquote('%E8%A1%97%E6%8B%8D'),
        'pd':'atlas',
        'dvpf':'pc',
        'aid':4916,
        'page_num':page_num,
        'search_json':'%7B%22from_search_id%22%3A%22202106100003510102121720341003A4ED%22%2C%22origin_keyword%22%3A%22%E8%A1%97%E6%8B%8D%22%2C%22image_keyword%22%3A%22%E8%A1%97%E6%8B%8D%22%7D',
        'rawJSON': 1,
        'search_id':'202106100004290101500200495C05B763'
    }
    url='https://so.toutiao.com/search?'+urlencode(params)
    try:
        response=requests.get(url,headers=headers,params=params)
        if response.status_code==200:
            return response.json()
    except requests.ConnectionError:
        return None

#获取图片链接
def get_images(json):
    images=json.get('rawData').get('data')
    for image in images:
        link = image.get('img_url')
        yield link

#下载图片
name=1
def saving_img(link):
    global name
    print(f'-------正在打印第{name}张图片')
    data=requests.get(link,headers=headers).content
    with open(f'image1/{name}.jpg','wb')as f:
            f.write(data)
            name+=1

def main(paga_num):
    json=get_page(paga_num)
    for link in get_images(json):
        saving_img(link)

if __name__ == '__main__':
    for i in range(0,2):
        main(i)