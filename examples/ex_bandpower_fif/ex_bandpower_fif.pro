#
#==============================================================================================================

include(../../mne-cpp.pri)

TEMPLATE = app

QT += widgets 3dextras

CONFIG   += console
CONFIG   -= app_bundle

TARGET = ex_bandpower_fif

CONFIG(debug, debug|release) {
    TARGET = $$join(TARGET,,,d)
}

DESTDIR =  $${MNE_BINARY_DIR}

contains(MNECPP_CONFIG, static) {
    CONFIG += static
    DEFINES += STATICBUILD
}

LIBS += -L$${MNE_LIBRARY_DIR}
CONFIG(debug, debug|release) {
    LIBS += -lmnecppRtProcessingd \
            -lmnecppMned \
            -lmnecppFiffd \
            -lmnecppUtilsd \
} else {
    LIBS +=-lmnecppRtProcessing \
            -lmnecppMne \
            -lmnecppFiff \
            -lmnecppUtils \
}

SOURCES += \
        main.cpp \

HEADERS += \

INCLUDEPATH += $${EIGEN_INCLUDE_DIR}
INCLUDEPATH += $${MNE_INCLUDE_DIR}

contains(MNECPP_CONFIG, withCodeCov) {
    QMAKE_CXXFLAGS += --coverage
    QMAKE_LFLAGS += --coverage
}

# Deploy dependencies
win32:!contains(MNECPP_CONFIG, static) {
    EXTRA_ARGS =
    DEPLOY_CMD = $$winDeployAppArgs($${TARGET},$${MNE_BINARY_DIR},$${MNE_LIBRARY_DIR},$${EXTRA_ARGS})
    QMAKE_POST_LINK += $${DEPLOY_CMD}
}
unix:!macx {
    # Unix
    QMAKE_RPATHDIR += $ORIGIN/../lib
}

# Activate FFTW backend in Eigen
contains(MNECPP_CONFIG, useFFTW):!contains(MNECPP_CONFIG, static) {
    DEFINES += EIGEN_FFTW_DEFAULT
    INCLUDEPATH += $$shell_path($${FFTW_DIR_INCLUDE})
    LIBS += -L$$shell_path($${FFTW_DIR_LIBS})

    win32 {
        # On Windows
        LIBS += -llibfftw3-3 \
                -llibfftw3f-3 \
                -llibfftw3l-3 \
    }

    unix:!macx {
        # On Linux
        LIBS += -lfftw3 \
                -lfftw3_threads \
    }
}
