from random import *        #导入随机数模板
from tkinter import *       #导入界面设计版块
import tkinter.font as tkFont #从tkinter.font模板导入字体设置函数


#自定义函数：运行后点击界面按钮进行的操作
def call1():
    
    r=randint(1,len(zhenxinhua)-2)      #randint（）函数表示随机取出一个范围在1到len(zhenxinhua)-2之间的随机数
    var.set(zhenxinhua[r])              #将取出的随机数作为序号把序号对应的题目放进var里

def call2():
    
    r=randint(1,len(damaoxian)-1)
    var.set(damaoxian[r])


#题库文件读取
try:
    f=open("真心话大冒险题库.txt")
    problem=[]                              #定义一个列表储存题库文件所有内容
    zhenxinhua=[]                           #定义列表储存真心话题目
    damaoxian=[]                            #定义列表储存大冒险题目
    for i in range(1000):                   #这个范围1000说明文件内容不能超过1000行，否则程序报错
        str=f.readline()                    #readline()函数可以从文件读取一行内容，并把其储存在str字符串里
        problem.append(str)                 #将str存入problem列表
        if(str==''):                        #碰到文件末尾时，直接结束循环
            break
    problem.remove('\n')                    #这两行是去掉读取到的无用内容
    problem.remove('')
    for i in range(1000):                   #这个循环是把problem里面真心话题目部分存进zhenxinhua列表
        if i>=len(problem):
            break
        zhenxinhua.append(problem[i])
        if problem[i]=='大冒险\n':
            break
for i in range(1000):                   #同上，这个是大冒险题目
    if problem[i]=='大冒险\n':
        while i :
            if i>=len(problem):
                break
                damaoxian.append(problem[i])
                i+=1
            break
except:
    
    print("error...")
    print("continue")



#界面创建主要操作
root = Tk()                                         #定义一个窗口，名字为root
root.title("真心话大冒险  v 2.01")                  #定义窗口显示出来时的标题
root.geometry('700x700')                            #用geometry固定窗口大小
root.attributes('-toolwindow',1)
ft=tkFont.Font(family='Comic Sans MS')              #定义一个字体变量
frame1=Frame(root,bg="green")                       #在root窗口里划分出一个区域frame1，定义背景颜色为绿色
frame2=Frame(root,bg="yellow")                      #在root窗口里划分出一个区域frame2，定义背景颜色为黄色
frame1.pack(padx=100,pady=100)                      #确定frame1的位置，参数表示padx表示区域与窗口竖边的距离，也表示与其他同类部件的距离，pady表示与横边的距离，也表示与其他同类部件的距离
frame2.pack(padx=30,pady=50)                               #确定frame2的位置，参数表示padx表示区域与窗口竖边的距离，也表示与其他同类部件的距离，pady表示与横边的距离，也表示与其他同类部件的距离
#上面隐藏了一个参数side，默认为TOP，表示居于窗口顶部中间


var=StringVar()                                     #一个字符串类，把字符串储存在var里，方便后面在窗口显示，显示的是题目
var.set("选择真心话或者大冒险")                     #设定初始字符串，set（）函数用于把参数写进var
textLabel=Label(frame1,textvariable=var,height=10, width=400,
                anchor=NW, wraplength=400,font=ft,bg="green")#Label，标签函数，用于在窗口中显示文字，将标签放在区域frame1里，调用var，将var里的字符串显示出来
textLabel.pack()                                    #确定标签的位置，同上面frame1.pack，只需要默认为居中

thebutton1=Button(frame2,text="真心话",width=15,height=1,font=ft,command=call1,bg="red")#在窗口设置一个按钮，命名为thebutton1，按钮上显示“真心话”，参数command用于当按钮被点击时调用的操作，bg表示按钮颜色
thebutton1.pack(side=LEFT,fill=Y,padx=50)           #确定按钮位置，放于frame2内


thebutton2=Button(frame2,text="大冒险",width=15,height=1,font=ft,command=call2,bg="blue")#同上
thebutton2.pack(side=RIGHT,fill=Y,padx=50)


mainloop()                                          #开始把创建好的窗口显示出来


