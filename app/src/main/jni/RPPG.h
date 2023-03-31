//
// Created by Mo3taz kayad on 2/9/2023.
//

#ifndef LIVENATIVERPPG_RPPG_H
#define LIVENATIVERPPG_RPPG_H


#include <fstream>
#include <string>
#include <stdio.h>
#include <opencv2/core.hpp>
#include <jni.h>

using namespace cv;
using namespace std;

enum RPPGAlgorithm { g, pca, xminay };

class RPPG {

public:

    // Constructor
    RPPG() {;}

    // Load Settings
    bool load(
//            jobject listener,
//              JNIEnv *jenv,                                   // Listener and environment for Java callback
              int algorithm,const double timeBase, const int downsample,
              const double samplingFrequency, const double rescanFrequency,
              const int minSignalSize, const int maxSignalSize,
              const string &logPath,
              const bool log, const bool gui);

    void processFrame(Mat &frameRGB, Mat &frameGray, int time);
    void processFrameWithNoFD(Mat &frameRGB , Mat &frameGray ,Mat& Mask, int time, double luma_val , double area);
    void assignArea(int &width, int &height);

    double getMeanBpm();
    void exit(JNIEnv *jenv);

    typedef vector<Point2f> Contour2f;

private:
    void extractSignal_g();
    void extractSignal_pca();
    void extractSignal_xminay();
    void estimateHeartrate(Mat &rgb);
    void invalidateFace();
    void log();

    void callback(int64_t now, double meanBpm, double minBpm, double maxBpm);   // Callback to Java

    void newPointGenerated(int64_t time, int p1 , int p2);

    // The JavaVM
//    JavaVM *jvm;

    // The listener
//    jobject listener;


    // The algorithm
    RPPGAlgorithm algorithm;

    // Settings
    Size minFaceSize;
    int maxSignalSize;
    int minSignalSize;
    double rescanFrequency;
    double samplingFrequency;
    double timeBase;
    bool logMode;
    bool guiMode;

    // State variables
    int64_t time;
    double fps;
    int high;
    int64_t lastSamplingTime;
    int64_t lastScanTime;
    int low;
    int64_t now;
    bool faceValid;
    bool rescanFlag;

    // Tracking
    Mat lastFrameGray;
    Contour2f corners;

    // Raw signal
    Mat1d s;
    Mat1d t;
    Mat1b re;

    // Estimation
    Mat1d s_f;
    Mat1d bpms;
    //Mat1d bpms_ws;
    Mat1d powerSpectrum;
    double bpm = 0.0;
    //double bpm_ws = 0.0;
    double meanBpm;
    double minBpm;
    double maxBpm;
    //double meanBpm_ws;
    //double minBpm_ws;
    //double maxBpm_ws;

    // Logfiles
    ofstream logfile;
    ofstream logfileDetailed;
    string logfilepath;


};



#endif //LIVENATIVERPPG_RPPG_H
