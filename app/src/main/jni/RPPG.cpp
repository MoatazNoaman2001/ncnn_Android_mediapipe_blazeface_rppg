//
// Created by Mo3taz kayad on 2/9/2023.
//

#include "RPPG.h"

#include <android/log.h>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/video/video.hpp>

#include "opencv.h"

using namespace cv;
using namespace std;

#define LOW_BPM 42
#define HIGH_BPM 240
#define REL_MIN_FACE_SIZE 0.2
#define SEC_PER_MIN 60
#define MAX_CORNERS 10
#define MIN_CORNERS 5
#define QUALITY_LEVEL 0.01
#define MIN_DISTANCE 25

#define LOG_TAG "Heartbeat::RPPG"
#define LOGD(...) ((void)__android_log_print(ANDROID_LOG_DEBUG, LOG_TAG, __VA_ARGS__))

bool RPPG::load(
//        jobject listener,
//        JNIEnv *jenv,
        int algorithm, const double timeBase, const int downsample,
        const double samplingFrequency, const double rescanFrequency,
        const int minSignalSize, const int maxSignalSize,
        const string &logPath,
        const bool log, const bool gui) {

    this->algorithm = (RPPGAlgorithm) algorithm;
    this->guiMode = gui;
    this->lastSamplingTime = 0;
    this->logMode = log;
    this->maxSignalSize = maxSignalSize;
    this->minSignalSize = minSignalSize;
    this->rescanFlag = false;
    this->rescanFrequency = rescanFrequency;
    this->samplingFrequency = samplingFrequency;
    this->timeBase = timeBase;

    LOGD("Using algorithm %d", algorithm);

    // Save reference to Java VM
//    jenv->GetJavaVM(&jvm);

    // Save global reference to listener object
//    this->listener = jenv->NewGlobalRef(listener);


    // Setting up logfilepath
    std::ostringstream path_1;
    path_1 << logPath << "_a=" << algorithm << "_min=" << minSignalSize << "_max=" << maxSignalSize
           << "_ds=" << downsample;
    this->logfilepath = path_1.str();

    // Logging bpm according to sampling frequency
    std::ostringstream path_2;
    path_2 << logfilepath << "_bpm.csv";
    logfile.open(path_2.str().c_str());
    logfile << "time;face_valid;mean;min;max\n";
    logfile.flush();

    // Logging bpm detailed
    std::ostringstream path_3;
    path_3 << logfilepath << "_bpmAll.csv";
    logfileDetailed.open(path_3.str().c_str());
    logfileDetailed << "time;face_valid;bpm\n";
    logfileDetailed.flush();

    return true;
}

void RPPG::exit(JNIEnv *jenv) {
//    jenv->DeleteGlobalRef(listener);
//    listener = NULL;
//    logfile.close();
//    logfileDetailed.close();
}

void RPPG::processFrameWithNoFD(Mat &frameRGB, Mat &frameGray, Mat &Mask, int time, double luma_val , double area) {
    this->time = time;


    if (luma_val > 150 && luma_val < 180 && area < 190 &&
        area > 100){
        lastScanTime = time;
        faceValid = true;
        rescanFlag = true;
    }else{
        lastScanTime = time;
//        invalidateFace();
    };

    if(true) {
        // Update fps
        fps = getFps(t, timeBase);
        LOGD("getFPS ok");

        // Remove old values from buffer
        while (s.rows > fps * maxSignalSize) {
            push(s);
            push(t);
            push(re);
        }
        LOGD("remove old values ok");

        assert(s.rows == t.rows && s.rows == re.rows);

        LOGD("assert 131 ok");

        // New values
        Scalar means = mean(frameRGB, Mask);

        LOGD("mean mask ok");
        // Add new values to raw signal buffer
        double values[] = {means(0), means(1), means(2)};
        s.push_back(Mat(1, 3, CV_64F, values));
        t.push_back(time);

        LOGD("add new values to raw signal buffer ok %d, %d , %d" , values[0],values[1],values[2]);

        // Save rescan flag
        re.push_back<bool>(rescanFlag);

        // Update fps
        fps = getFps(t, timeBase);
        // Update band spectrum limits
        LOGD("s rows: %i , fps: %f , timeBase: %f", s.rows, fps, timeBase);
        low = (int) (s.rows * LOW_BPM / SEC_PER_MIN / fps);
        high = (int) (s.rows * HIGH_BPM / SEC_PER_MIN / fps) + 1;
        LOGD("159 ok low: %i , high:%i ", low, high);


        // If valid signal is large enough: estimate
        if (s.rows >= fps * minSignalSize) {
            // Filtering
            switch (algorithm) {
                case g:
                    extractSignal_g();
                    break;
                case pca:
                    extractSignal_pca();
                    break;
                case xminay:
                    extractSignal_xminay();
                    break;
            }

            // PSD estimation
            estimateHeartrate(frameRGB);

            // Log
            log();
        }
    }

    rescanFlag = false;

    frameGray.copyTo(lastFrameGray);
}

void RPPG::processFrame(Mat &frameRGB, Mat &frameGray, int time) {

//    // Set time
//    this->time = time;
//
//    if (!faceValid) {
//
//        LOGD("Not valid, finding a new face");
//
//        lastScanTime = time;
//        detectFace(frameRGB, frameGray);
//
//    } else if ((time - lastScanTime) * timeBase >= 1 / rescanFrequency) {
//
//        LOGD("Valid, but rescanning face");
//
//        lastScanTime = time;
//        detectFace(frameRGB, frameGray);
//        rescanFlag = true;
//
//    } else {
//
//        LOGD("Tracking face");
//
//        trackFace(frameGray);
//    }
//
//    if (faceValid) {
//
//        // Update fps
//        fps = getFps(t, timeBase);
//        LOGD("getFPS ok");
//
//        // Remove old values from buffer
//        while (s.rows > fps * maxSignalSize) {
//            push(s);
//            push(t);
//            push(re);
//        }
//        LOGD("remove old values ok");
//
//        assert(s.rows == t.rows && s.rows == re.rows);
//
//        LOGD("assert 131 ok");
//
//        // New values
//        Scalar means = mean(frameRGB, mask);
//
//        LOGD("mean mask ok");
//        // Add new values to raw signal buffer
//        double values[] = {means(0), means(1), means(2)};
//        s.push_back(Mat(1, 3, CV_64F, values));
//        t.push_back(time);
//
//        LOGD("add new values to raw signal buffer ok");
//
//        // Save rescan flag
//        re.push_back<bool>(rescanFlag);
//
//        // Update fps
//        fps = getFps(t, timeBase);
//        // Update band spectrum limits
//        LOGD("s rows: %i , fps: %f , timeBase: %f" , s.rows , fps , timeBase);
//        low = (int) (s.rows * LOW_BPM / SEC_PER_MIN / fps);
//        high = (int) (s.rows * HIGH_BPM / SEC_PER_MIN / fps) + 1;
//        LOGD("159 ok low: %i , high:%i " , low, high);
//
//
//        // If valid signal is large enough: estimate
//        if (s.rows >= fps * minSignalSize) {
//            // Filtering
//            switch (algorithm) {
//                case g:
//                    extractSignal_g();
//                    break;
//                case pca:
//                    extractSignal_pca();
//                    break;
//                case xminay:
//                    extractSignal_xminay();
//                    break;
//            }
//
//            // PSD estimation
//            estimateHeartrate();
//
//            // Log
//            log();
//        }
//
//        if (guiMode) {
//            draw(frameRGB);
//        }
//    }
//
//    if (!guiMode) {
//        // Indicator
//        frameRGB.setTo(BLACK);
//        // circle(frameRGB, Point(1250, 100), 25, faceValid ? GREEN : RED, -1, 8, 0);
//    }
//
//    rescanFlag = false;
//
//    frameGray.copyTo(lastFrameGray);
}


void RPPG::invalidateFace() {

    s = Mat1d();
    s_f = Mat1d();
    t = Mat1d();
    re = Mat1b();
    powerSpectrum = Mat1d();
    faceValid = false;
}

void RPPG::extractSignal_g() {

    // Denoise
    Mat s_den = Mat(s.rows, 1, CV_64F);
    denoise(s.col(1), re, s_den);
    LOGD("Denoise ok");

    // Normalise
    normalization(s_den, s_den);
    LOGD("Normalise ok");

    // Detrend
    Mat s_det = Mat(s_den.rows, s_den.cols, CV_64F);
    detrend(s_den, s_det, fps);
    LOGD("Detrend ok");

    // Moving average
    Mat s_mav = Mat(s_det.rows, s_det.cols, CV_64F);
    movingAverage(s_det, s_mav, 3, fmax(floor(fps / 6), 2));
    LOGD("moving avg ok");

    s_mav.copyTo(s_f);

    // Logging
    if (logMode) {
        std::ofstream log;
        std::ostringstream filepath;
        filepath << logfilepath << "_signal_" << time << ".csv";
        log.open(filepath.str().c_str());
        log << "g;g_den;g_det;g_mav\n";
        for (int i = 0; i < s.rows; i++) {
            log << s.at<double>(i, 1) << ";";
            log << s_den.at<double>(i, 0) << ";";
            log << s_det.at<double>(i, 0) << ";";
            log << s_mav.at<double>(i, 0) << "\n";
        }
        log.close();
    }
    LOGD("logging ok");
}

void RPPG::extractSignal_pca() {

    // Denoise signals
    Mat s_den = Mat(s.rows, s.cols, CV_64F);
    denoise(s, re, s_den);

    // Normalize signals
    normalization(s_den, s_den);

    // Detrend
    Mat s_det = Mat(s.rows, s.cols, CV_64F);
    detrend(s_den, s_det, fps);

    // PCA to reduce dimensionality
    Mat s_pca = Mat(s.rows, 1, CV_32F);
    Mat pc = Mat(s.rows, s.cols, CV_32F);
    pcaComponent(s_det, s_pca, pc, low, high);

    // Moving average
    Mat s_mav = Mat(s.rows, 1, CV_32F);
    movingAverage(s_pca, s_mav, 3, fmax(floor(fps / 6), 2));

    s_mav.copyTo(s_f);

    // Logging
    if (logMode) {
        std::ofstream log;
        std::ostringstream filepath;
        filepath << logfilepath << "_signal_" << time << ".csv";
        log.open(filepath.str().c_str());
        log << "re;r;g;b;r_den;g_den;b_den;r_det;g_det;b_det;pc1;pc2;pc3;s_pca;s_mav\n";
        for (int i = 0; i < s.rows; i++) {
            log << re.at<bool>(i, 0) << ";";
            log << s.at<double>(i, 0) << ";";
            log << s.at<double>(i, 1) << ";";
            log << s.at<double>(i, 2) << ";";
            log << s_den.at<double>(i, 0) << ";";
            log << s_den.at<double>(i, 1) << ";";
            log << s_den.at<double>(i, 2) << ";";
            log << s_det.at<double>(i, 0) << ";";
            log << s_det.at<double>(i, 1) << ";";
            log << s_det.at<double>(i, 2) << ";";
            log << pc.at<double>(i, 0) << ";";
            log << pc.at<double>(i, 1) << ";";
            log << pc.at<double>(i, 2) << ";";
            log << s_pca.at<double>(i, 0) << ";";
            log << s_mav.at<double>(i, 0) << "\n";
        }
        log.close();
    }
}

void RPPG::extractSignal_xminay() {

    // Denoise signals
    Mat s_den = Mat(s.rows, s.cols, CV_64F);
    denoise(s, re, s_den);

    // Normalize raw signals
    Mat s_n = Mat(s_den.rows, s_den.cols, CV_64F);
    normalization(s_den, s_n);

    // Calculate X_s signal
    Mat x_s = Mat(s.rows, s.cols, CV_64F);
    addWeighted(s_n.col(0), 3, s_n.col(1), -2, 0, x_s);

    // Calculate Y_s signal
    Mat y_s = Mat(s.rows, s.cols, CV_64F);
    addWeighted(s_n.col(0), 1.5, s_n.col(1), 1, 0, y_s);
    addWeighted(y_s, 1, s_n.col(2), -1.5, 0, y_s);

    // Bandpass
    Mat x_f = Mat(s.rows, s.cols, CV_32F);
    bandpass(x_s, x_f, low, high);
    x_f.convertTo(x_f, CV_64F);
    Mat y_f = Mat(s.rows, s.cols, CV_32F);
    bandpass(y_s, y_f, low, high);
    y_f.convertTo(y_f, CV_64F);

    // Calculate alpha
    Scalar mean_x_f;
    Scalar stddev_x_f;
    meanStdDev(x_f, mean_x_f, stddev_x_f);
    Scalar mean_y_f;
    Scalar stddev_y_f;
    meanStdDev(y_f, mean_y_f, stddev_y_f);
    double alpha = stddev_x_f.val[0] / stddev_y_f.val[0];

    // Calculate signal
    Mat xminay = Mat(s.rows, 1, CV_64F);
    addWeighted(x_f, 1, y_f, -alpha, 0, xminay);

    // Moving average
    movingAverage(xminay, s_f, 3, fmax(floor(fps / 6), 2));

    // Logging
    if (logMode) {
        std::ofstream log;
        std::ostringstream filepath;
        filepath << logfilepath << "_signal_" << time << ".csv";
        log.open(filepath.str().c_str());
        log << "r;g;b;r_den;g_den;b_den;x_s;y_s;x_f;y_f;s;s_f\n";
        for (int i = 0; i < s.rows; i++) {
            log << s.at<double>(i, 0) << ";";
            log << s.at<double>(i, 1) << ";";
            log << s.at<double>(i, 2) << ";";
            log << s_den.at<double>(i, 0) << ";";
            log << s_den.at<double>(i, 1) << ";";
            log << s_den.at<double>(i, 2) << ";";
            log << x_s.at<double>(i, 0) << ";";
            log << y_s.at<double>(i, 0) << ";";
            log << x_f.at<double>(i, 0) << ";";
            log << y_f.at<double>(i, 0) << ";";
            log << xminay.at<double>(i, 0) << ";";
            log << s_f.at<double>(i, 0) << "\n";
        }
        log.close();
    }
}

void RPPG::estimateHeartrate(Mat &frameRGB) {
    LOGD("estimate Heart rate bigin");
    powerSpectrum = Mat(s_f.size(), s_f.type());
    timeToFrequency(s_f, powerSpectrum, true);
    LOGD("time to freq ok");

    // band mask
    const int total = s_f.rows;
    Mat bandMask = Mat::zeros(s_f.size(), CV_8U);
    bandMask.rowRange(min(low, total), min(high, total) + 1).setTo(ONE);


    LOGD("band mask ok");


    if (!powerSpectrum.empty()) {

        // grab index of max power spectrum
        double min, max;
        Point pmin, pmax;
        minMaxLoc(powerSpectrum, &min, &max, &pmin, &pmax, bandMask);

        // calculate BPM
        bpm = pmax.y * fps / total * SEC_PER_MIN;
        bpms.push_back(bpm);

        // calculate BPM based on weighted squares power spectrum
        //double weightedSquares = weightedSquaresMeanIndex(powerSpectrum, low, high);
        //double bpm_ws = weightedSquares * fps / total * SEC_PER_MIN;
        //bpms_ws.push_back(bpm_ws);

        LOGD("FPS=%f Vals=%d Peak=%d BPM=%f", fps, powerSpectrum.rows, pmax.y, bpm);

        // Logging
        if (logMode) {
            std::ofstream log;
            std::ostringstream filepath;
            filepath << logfilepath << "_estimation_" << time << ".csv";
            log.open(filepath.str().c_str());
            log << "i;powerSpectrum\n";
            for (int i = 0; i < powerSpectrum.rows; i++) {
                if (low <= i && i <= high) {
                    log << i << ";";
                    log << powerSpectrum.at<double>(i, 0) << "\n";
                }
            }
            log.close();
        }
    }

    if ((time - lastSamplingTime) * timeBase >= 1 / samplingFrequency) {
        lastSamplingTime = time;

        cv::sort(bpms, bpms, SORT_EVERY_COLUMN);

        // average calculated BPMs since last sampling time
        meanBpm = mean(bpms)(0);
        minBpm = bpms.at<double>(0, 0);
        maxBpm = bpms.at<double>(bpms.rows - 1, 0);

        // cv::sort(bpms_ws, bpms_ws, SORT_EVERY_COLUMN);
        // meanBpm_ws = mean(bpms_ws)(0);
        // minBpm_ws = bpms_ws.at<double>(0, 0);
        // maxBpm_ws = bpms_ws.at<double>(bpms_ws.rows-1, 0);

        callback(time, meanBpm, minBpm, maxBpm);

        bpms.pop_back(bpms.rows);
        // bpms_ws.pop_back(bpms_ws.rows);
    }
}

void RPPG::log() {

    if (lastSamplingTime == time || lastSamplingTime == 0) {
        logfile << time << ";";
        logfile << faceValid << ";";
        logfile << meanBpm << ";";
        logfile << minBpm << ";";
        logfile << maxBpm << "\n";
        logfile.flush();
    }

    logfileDetailed << time << ";";
    logfileDetailed << faceValid << ";";
    logfileDetailed << bpm << "\n";
    logfileDetailed.flush();
}

void RPPG::newPointGenerated(int64_t time, int p1, int p2) {

//    JNIEnv *jenv;
//    int stat = jvm->GetEnv((void **) &jenv, JNI_VERSION_1_6);
//
//    if (stat == JNI_EDETACHED) {
//        LOGD("GetEnv: not attached");
//        if (jvm->AttachCurrentThread(&jenv, NULL) != 0) {
//            LOGD("GetEnv: Failed to attach");
//        } else {
//            LOGD("GetEnv: Attached to %i", jenv);
//        }
//    } else if (stat == JNI_OK) {
//        //
//    } else if (stat == JNI_EVERSION) {
//        LOGD("GetEnv: version not supported");
//    }

//    jclass classRef = jenv->FindClass("com/example/livenativerppg/component/natives/SignalPoint");
//    jmethodID mid = jenv->GetMethodID(classRef , "<init>", "(JII)V");
//    jobject cons_obj = jenv->NewObject(classRef , mid , time , p1, p2);

//    jclass listenerClassRef = jenv->GetObjectClass(listener);
//    jmethodID methodId = jenv->GetMethodID(listenerClassRef , "onNewPointGenerated" , "(Lcom/example/livenativerppg/component/natives/SignalPoint;)V");
//    jenv->CallVoidMethod(listenerClassRef , methodId , cons_obj);

//    jenv->DeleteLocalRef(cons_obj);
}

void RPPG::callback(int64_t time, double meanBpm, double minBpm, double maxBpm) {

//    JNIEnv *jenv;
//    int stat = jvm->GetEnv((void **) &jenv, JNI_VERSION_1_6);
//
//    if (stat == JNI_EDETACHED) {
//        LOGD("GetEnv: not attached");
//        if (jvm->AttachCurrentThread(&jenv, NULL) != 0) {
//            LOGD("GetEnv: Failed to attach");
//        } else {
//            LOGD("GetEnv: Attached to %d", jenv);
//        }
//    } else if (stat == JNI_OK) {
//        //
//    } else if (stat == JNI_EVERSION) {
//        LOGD("GetEnv: version not supported");
//    }

//    // Return object
//
//    // Get Return object class reference
//    jclass returnObjectClassRef = jenv->FindClass(
//            "com/example/livenativerppg/component/natives/RPPGResult");
//
//    // Get Return object constructor method
//    jmethodID constructorMethodID = jenv->GetMethodID(returnObjectClassRef, "<init>", "(JDDD)V");
//
//    // Create Info class
//    jobject returnObject = jenv->NewObject(returnObjectClassRef, constructorMethodID, (jlong) time,
//                                           meanBpm, minBpm, maxBpm);
//
//    // Listener
//
//    // Get the Listener class reference
//    jclass listenerClassRef = jenv->GetObjectClass(listener);
//
//    // Use Listener class reference to load the eventOccurred method
//    jmethodID listenerEventOccuredMethodID = jenv->GetMethodID(listenerClassRef, "onRPPGResult", "(Lcom/example/livenativerppg/component/natives/RPPGResult;)V");
//
//    // Invoke listener eventOccurred
//    jenv->CallVoidMethod(listener, listenerEventOccuredMethodID, returnObject);
//
//    // Cleanup
//    jenv->DeleteLocalRef(returnObject);
}


void RPPG::assignArea(int &width, int &height) {

}

double RPPG::getMeanBpm() {
    return meanBpm;
}

Mat1d RPPG::s_f_return() {
    return s_f;
}
