#include <iostream>
#include <direct.h>
#include <cstdio>
#include <iomanip>
#include <string>
#include "opencv2\opencv.hpp"
#include <fstream>

using namespace std;
using namespace cv;

string filename="stainedcell_small";	//檔名
string input_type=".jpg";				//檔案類型	
int Threshold=190;						//二值化閾值，small時用195，large時用190
int Scalar_slider=1;					//放大倍率
int closest;							//最近點編號
Rect bounding;							//目標細胞範圍

const int Scalar_slider_max = 10;		//tracker的最大值
const double pi=3.14159;

Mat frame,frame_wshd,frame_target,target;	//frame為彩圖 frame_wshd存watershed結果 frame_target為目標細胞 target為目標細胞放大後
Point point_mouse;					//紀錄目前滑鼠游標點

class WatershedSegment{
private:
    Mat markers;
public:
    void setMarkers(Mat &markerImage)
    {
        markerImage.convertTo(markers,CV_32S);//convert to image of ints
 
    }
 
    Mat process(Mat &image)
    {
        watershed(image,markers);
        markers.convertTo(markers,CV_8U);
        return markers;
    }
};

Mat wshd(Mat img){
	Mat img_gray,img_fore,img_back,result;
	Mat markers(img_gray.size(),CV_8U,Scalar(0));
	WatershedSegment segmenter;
	vector<vector<Point> > contour;
    vector<Vec4i> hierarchy;
	
	cvtColor(img, img_gray, COLOR_BGR2GRAY);
	threshold(img_gray,img_gray,Threshold,255,THRESH_BINARY);
	
	//foreground
	erode(img_gray,img_fore,Mat(),Point(-1,-1),2);

	//background
	dilate(img_gray,img_back,Mat(),Point(-1,-1),3);
    threshold(img_back,img_back,1,128,THRESH_BINARY_INV);

	//creat markers
	markers=img_fore+img_back;

	//watershed
	segmenter.setMarkers(markers);
	result=segmenter.process(img);
	threshold(result,result,Threshold,255,THRESH_BINARY);

	return result;
}

void discription(Mat result,Mat img,bool initial=true){
	int N;
	double area_tot=0,area;
	double len_tot=0,len;
    double compactness=0,regularity=0;
	double var,d,d_avg;
	string comp,reg;
	Point point,point_closest;
	vector<Vec4i> hierarchy;
	static vector<vector<Point> > contour;
	
	if(initial){
		threshold(result,result,Threshold,255,THRESH_BINARY);
		findContours(result, contour, hierarchy, RETR_TREE, CV_CHAIN_APPROX_SIMPLE,Point(0,0));	//search contour
	}
	vector<Moments> mu(contour.size());	//moments
	static vector<Point2f> mc(contour.size());	//centers
	
	if(initial==true){
		
		ofstream outputfile(filename+"//"+filename+"_discription_result"+".txt");

		if( contour.empty() )
			cout<<"contour is empty!!"<<endl;
		else{
			for(int i=1;i<contour.size();i++){
				var=0;
				d_avg=0;
				area=contourArea(contour[i],false);		//calculate area
				len=arcLength(contour[i],true);			//calculate closed curve length
				mu[i] = moments( contour[i], false );	// Get the moments
			
				// Get the mass centers
				mc[i] = Point2f( mu[i].m10/mu[i].m00 , mu[i].m01/mu[i].m00 );	
				point.x=(int)mc[i].x;
				point.y=(int)mc[i].y;
			
				//calculate compactness
				compactness=(len*len)/(4*pi*area);		
				comp=(compactness>4)?"不算集中":(compactness>2)?"還算集中":"非常集中";
			
				//calculate regularity
				for(int j=0;j<contour[i].size();j++){
					d=(contour[i][j].x-mc[i].x)*(contour[i][j].x-mc[i].x)+(contour[i][j].y-mc[i].y)*(contour[i][j].y-mc[i].y);
					d_avg+=sqrt(d);
					var+=d;	
				}
				N=contour[i].size();
				d_avg/=N;
				regularity=sqrt((var-d_avg*d_avg*N)/N)/d_avg;
				reg=(regularity>0.2)?"不算規律":(regularity>0.1)?"還算規律":"非常規律";
			
				//output
				cout<<"#"<<right<<setw(3) <<i<<" Area = "<<left<<setw(6) <<area<<" pixels, "<<" Length = "<<left<<setw(7) <<len<< " pixels, "
					<<"Center at ("<<left<<setw(7)<<mc[i].x<<" , "<<left<<setw(7)<<mc[i].y<<" ), "
					<<"Compactness ="<<left<<setw(7)<<compactness<<" ("<<comp<<")"
					<<", Regularity ="<<left<<setw(9)<<regularity<<" ("<<reg<<")"<<endl;
				outputfile<<"#"<<right<<setw(3) <<i<<" Area = "<<left<<setw(6) <<area<<" pixels, "<<" Length = "<<left<<setw(7) <<len<< " pixels, "
					<<"Center at ("<<left<<setw(7)<<mc[i].x<<" , "<<left<<setw(7)<<mc[i].y<<" ), "
					<<"Compactness ="<<left<<setw(7)<<compactness<<" ("<<comp<<")"
					<<", Regularity ="<<left<<setw(9)<<regularity<<" ("<<reg<<")"<<endl;
				drawContours(img, contour,i, Scalar(0, 255, 0), 1);		//circule the cell
				circle(img, point, 2, Scalar(1, 247, 255), -1);			//circule the cell center
				putText(img,to_string(i), point, FONT_HERSHEY_SIMPLEX, 0.4, Scalar(0, 0, 255), 1);		//label the cell
			
				area_tot+=area;
				len_tot+=len;
			}
			outputfile.close();
			cout<<endl;
			cout<<"The average area is "<<area_tot/contour.size()<< " pixels"<<endl;
			cout<<"The average length is "<<len_tot/contour.size()<< " pixels"<<endl;
		}
		outputfile.close();
		imshow("Finish",img);
		imwrite(filename+"//"+filename+"_discription"+input_type,img);
	}
	if(initial==false){
		Mat showpic;
		frame.copyTo(showpic);
		
		static string str1,str2,str3,str4,str5,str6;
		double min=img.rows*img.rows+img.cols*img.cols;
		string Area,Len,Compactness,Regularity;
		
		closest=1;
		for(int i=1;i<contour.size();i++){
			d =(point_mouse.x-mc[i].x)*(point_mouse.x-mc[i].x)+(point_mouse.y-mc[i].y)*(point_mouse.y-mc[i].y);	
			if(d<min){
				closest=i;
				min=d;
			}
		}
		bounding=boundingRect(contour[closest]);
		point_closest.x=(int)mc[closest].x;
		point_closest.y=(int)mc[closest].y;
		area=contourArea(contour[closest],false);		//calculate area
		len=arcLength(contour[closest],true);			//calculate closed curve length
		
		//calculate compactness
		compactness=(len*len)/(4*pi*area);		
		comp=(compactness>4)?"Bad":(compactness>2)?"Fair":"Good";
		
		//calculate regularity
		var=0;
		d_avg=0;
		for(int j=0;j<contour[closest].size();j++){
			d=(contour[closest][j].x-mc[closest].x)*(contour[closest][j].x-mc[closest].x)+(contour[closest][j].y-mc[closest].y)*(contour[closest][j].y-mc[closest].y);
			d_avg+=sqrt(d);
			var+=d;	
		}
		N=contour[closest].size();
		d_avg/=N;
		regularity=sqrt((var-d_avg*d_avg*N)/N)/d_avg;
		reg=(regularity>0.2)?"Bad":(regularity>0.1)?"Fair":"Good";
		
		Area.assign(to_string(area),0,7);
		Len.assign(to_string(len),0,7);
		Compactness.assign(to_string(compactness),0,7);
		Regularity.assign(to_string(regularity),0,7);

		str1="#"+to_string(closest);
		str2="Area = "+Area+" pixels";
		str3="Length = "+Len+ " pixels";
		str4="Center at ("+to_string(mc[closest].x)+" , "+to_string(mc[closest].y)+" )";
		str5="Compactness ="+Compactness+" ("+comp+")";
		str6="Regularity ="+Regularity+" ("+reg+")";

		drawContours(showpic, contour,closest, Scalar(0, 255, 0), 1);		//circule the cell
		circle(showpic, point_closest, 2, Scalar(1, 247, 255), -1);			//circule the cell center
		putText(showpic,to_string(closest), point_closest, FONT_HERSHEY_SIMPLEX, 0.4, Scalar(0, 0, 255), 1);		//label the cell
		putText(showpic,str1,Point(10,10), FONT_HERSHEY_SIMPLEX, 0.3, Scalar(255, 0, 0), 1);		//show the information of the target cell
		putText(showpic,str2,Point(10,20), FONT_HERSHEY_SIMPLEX, 0.3, Scalar(255, 0, 0), 1);		//show the information of the target cell
		putText(showpic,str3,Point(10,30), FONT_HERSHEY_SIMPLEX, 0.3, Scalar(255, 0, 0), 1);		//show the information of the target cell
		putText(showpic,str4,Point(10,40), FONT_HERSHEY_SIMPLEX, 0.3, Scalar(255, 0, 0), 1);		//show the information of the target cell
		putText(showpic,str5,Point(10,50), FONT_HERSHEY_SIMPLEX, 0.3, Scalar(255, 0, 0), 1);		//show the information of the target cell
		putText(showpic,str6,Point(10,60), FONT_HERSHEY_SIMPLEX, 0.3, Scalar(255, 0, 0), 1);		//show the information of the target cell

		imshow("original",showpic);
		
		frame(Rect(bounding.x,bounding.y,bounding.width, bounding.height)).copyTo(frame_target);
		
		Size dsize=Size(int(Scalar_slider*frame_target.cols), int(Scalar_slider*frame_target.rows));
		resize(frame_target,target,dsize,Scalar_slider,Scalar_slider,CV_INTER_CUBIC);
		imshow("target",target);
	}
}

void onMouse(int Event,int x,int y,int flags,void* param){			//滑鼠點擊事件
    if(Event==CV_EVENT_LBUTTONDOWN){		//右鍵為選取目標細胞
		//紀錄點擊座標
		point_mouse.x=x;
		point_mouse.y=y;
		
		destroyWindow("target");
		namedWindow("target",WINDOW_AUTOSIZE);
		
		discription(frame_wshd,frame,false);	//Discription

		}
	else if(Event==CV_EVENT_RBUTTONDOWN){		//左鍵為存檔功能
		imwrite(filename+"//"+filename+"_target_"+to_string(closest)+ input_type,frame_target);
	}
	else
		waitKey(1);
}

void on_trackbar( int, void* )
{
	int scalar=(Scalar_slider>1)?Scalar_slider:1;
	

	//把target的顯示視窗歸零
	destroyWindow("target");
	namedWindow("target",WINDOW_AUTOSIZE);

	//調整輸出視窗與目標細胞尺寸
	Size dsize=Size(int(scalar*frame_target.cols), int(scalar*frame_target.rows));
	resize(frame_target,target,dsize,scalar,scalar,CV_INTER_CUBIC);
	resizeWindow("target",dsize.width,dsize.height);
	
	//輸出目標細胞
	imshow("target",target);
}

void main(){
	Mat pic;		//用來輸出顯示的
	char folderpath[20];	//輸出存檔資料夾名

	frame=imread(filename+input_type,CV_LOAD_IMAGE_COLOR);
	frame.copyTo(pic);
	namedWindow("original",WINDOW_NORMAL);
	namedWindow("Finish",WINDOW_NORMAL);
	namedWindow("target",WINDOW_AUTOSIZE);
	imshow("original",frame);

	//創資料夾
	strcpy_s(folderpath,filename.c_str());	
	_mkdir(folderpath);
		
	//設定滑鼠點擊事件
	setMouseCallback("original",onMouse,NULL);
	
	createTrackbar("scalar","original",&Scalar_slider,Scalar_slider_max,on_trackbar);

	frame_wshd=wshd(frame);			//Watershed segmentation
	discription(frame_wshd,pic);	//Discription
	
	waitKey(0);
	system("pause");
}
