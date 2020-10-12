#include <iostream>
#include <fstream>
#include <string>
using namespace std;

Int_t makeLatexTextForPhotos_v(){

  //////////////////////////////////////////////////
  TString photoFolder = "photos_v";
  TString photoList   = "listOfPhotos.dat";
  TString outputFile  = "photos_v.tex";
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
