# OpenCV
win32{
    INCLUDEPATH *= C:/Development/opencv/build/include
    LIBS *= -L"C:/Development/opencv/build/x64/vc11/lib"
    OpenCV_VERSION = 249
    Debug:LIBS *= -lopencv_core$${OpenCV_VERSION}d -lopencv_highgui$${OpenCV_VERSION}d -lopencv_features2d$${OpenCV_VERSION}d -lopencv_nonfree$${OpenCV_VERSION}d -lopencv_flann$${OpenCV_VERSION}d -lopencv_imgproc$${OpenCV_VERSION}d
    Release:LIBS *= -lopencv_core$${OpenCV_VERSION} -lopencv_highgui$${OpenCV_VERSION} -lopencv_features2d$${OpenCV_VERSION} -lopencv_nonfree$${OpenCV_VERSION} -lopencv_flann$${OpenCV_VERSION} -lopencv_imgproc$${OpenCV_VERSION}
}
