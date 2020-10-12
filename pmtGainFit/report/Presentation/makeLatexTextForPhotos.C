#include <iostream>
#include <fstream>
#include <string>
using namespace std;

TString getDate(TString mot);
TString findStrAndReplace(string mot,string replOf,string replWith);

Int_t makeLatexTextForPhotos(){

  //////////////////////////////////////////////////
  TString photoFolder = "photos";
  TString photoList   = "listOfPhotos.dat";
  TString outputFile  = "photos.tex";
  Int_t nMaxPhotos    = 1000;
  //////////////////////////////////////////////////

  TString photoName[nMaxPhotos];
  TString photoTime[nMaxPhotos];
  TString photoList_file = photoFolder; photoList_file += "/"; photoList_file += photoList;
  TString mot;
  //cout<<"photoList_file "<<photoList_file<<endl;
  ifstream myfilePhotoList (photoList_file.Data());
  Int_t nPhotos = 0;
  Int_t i = 0;

  cout<<"photoFolder    = "<<photoFolder<<endl
      <<"photoList      = "<<photoList<<endl
      <<"photoList_file = "<<photoList_file<<endl
      <<"outputFile     = "<<outputFile<<endl;
  if(myfilePhotoList.is_open()){
    while ( myfilePhotoList>>mot ){
      //cout << mot << '\n';
      if(nPhotos<nMaxPhotos){
	photoTime[nPhotos] = getDate(mot);
	photoName[nPhotos] = "./"; 
	photoName[nPhotos] += photoFolder; 
	photoName[nPhotos] += "/"; 
	photoName[nPhotos] += mot; 
      }
      else{
	cout<<" ERROR --> nPhotos >= nMaxPhotos "<<endl
	    <<"           nPhotos  = "<<nPhotos<<endl
	    <<"        nMaxPhotos  = "<<nMaxPhotos<<endl;
	return -1;
      }
      nPhotos++;
    }
    myfilePhotoList.close();
  }
  else{
    cout << "Unable to open file"; 
  }
  //cout<<"nPhotos = "<<nPhotos<<endl;

  ofstream myfileOut;
  myfileOut.open(outputFile.Data());
  for(i = 0;i<nPhotos;i++){
    //%------------------------------------------------
    //\begin{frame}
    //\frametitle{Experiment description}
    //$2015.04.13-17:39:40$ 
    //\begin{figure}
    //\includegraphics[width=0.95\linewidth]{./photos/2015-04-13_17-39-40.jpg}
    //\end{figure}
    //\end{frame}
    //%------------------------------------------------
    myfileOut<<"%------------------------------------------------"<<endl;
    myfileOut<<"\\begin{frame}"<<endl;
    myfileOut<<"\\frametitle{Experiment description}"<<endl;
    myfileOut<<"$"<<photoTime[i]<<"$"<<endl;
    myfileOut<<"\\includegraphics[width=0.95\\linewidth]{"<<photoName[i]<<"}"<<endl;
    myfileOut<<"\\begin{figure}"<<endl;
    myfileOut<<"\\end{figure}"<<endl;
    myfileOut<<"\\end{frame}"<<endl;
    myfileOut<<"%------------------------------------------------"<<endl;
    myfileOut<<endl<<endl;
    //cout<<photoName[i]<<" "<<photoTime[i]<<endl;
  }
  myfileOut.close();

  return 0;
}

TString getDate(TString mot){
  TString motDate;
  string s = mot;
  string date_s;
  string time_s;
  TString date_ts;
  TString time_ts;
  size_t f = s.find(".jpg");
  //cout<<s<<endl;
  s.replace(f,std::string(".jpg").length(),"");
  //cout<<s<<endl;
  size_t f = s.find("_");
  s.replace(f,std::string("_").length()," ");
  //cout<<s<<endl;
  size_t f = s.find(" ");
  date_s = s.substr(0,f);
  time_s = s.substr(f+1,s.length());
  date_ts = findStrAndReplace(date_s,"-",".");
  time_ts = findStrAndReplace(time_s,"-",":");
  //cout<<date_s<<endl;
  //cout<<time_s<<endl;
  //cout<<date_ts<<endl;
  //cout<<time_ts<<endl;
  motDate = date_ts;
  motDate += "-";
  motDate += time_ts;
  //cout<<motDate<<endl;
  return motDate;
}

TString findStrAndReplace(string mot,string replOf,string replWith){
  string s = mot;
  TString dateTime;
  size_t f = s.find(replOf.data());
  s.replace(f,std::string(replOf.data()).length(),replWith.data());
  size_t f = s.find(replOf.data());
  s.replace(f,std::string(replOf.data()).length(),replWith.data());
  dateTime = s;
  return dateTime;
}
