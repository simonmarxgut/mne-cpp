//=============================================================================================================
/**
 * @file     main.cpp
 * @author   Johannes Vorwerk <johannes.vorwerk@umit.at>
 * @since    0.1.0
 * @date     11, 2019
 *
 * @section  LICENSE
 *
 * Copyright (C) 2021, Johannes Vorwerk. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that
 * the following conditions are met:
 *     * Redistributions of source code must retain the above copyright notice, this list of conditions and the
 *       following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
 *       the following disclaimer in the documentation and/or other materials provided with the distribution.
 *     * Neither the name of MNE-CPP authors nor the names of its contributors may be used
 *       to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * @brief     Example for filtering data with a user defined FIR filter and writing the result to a file.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include <iostream>
#include <vector>
#include <math.h>

#include <fiff/fiff.h>

#include <utils/filterTools/filterdata.h>
#include <utils/generics/applicationlogger.h>

#include <rtprocessing/rtfilter.h>

#include <utils/spectral.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QtCore/QCoreApplication>
#include <QFile>
#include <QCommandLineParser>

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace FIFFLIB;
using namespace UTILSLIB;
using namespace RTPROCESSINGLIB;
using namespace Eigen;

//=============================================================================================================
// MAIN
//=============================================================================================================

//=============================================================================================================
/**
 * The function main marks the entry point of the program.
 * By default, main has the storage class extern.
 *
 * @param [in] argc (argument count) is an integer that indicates how many arguments were entered on the command line when the program was started.
 * @param [in] argv (argument vector) is an array of pointers to arrays of character objects. The array objects are null-terminated strings, representing the arguments that were entered on the command line when the program was started.
 * @return the value that was set to exit() (which is 0 if exit() is called via quit()).
 */
int main(int argc, char *argv[])
{
    qInstallMessageHandler(ApplicationLogger::customLogWriter);
    QCoreApplication a(argc, argv);

    // Command Line Parser
    QCommandLineParser parser;
    parser.setApplicationDescription("Read Write Raw Example");
    parser.addHelpOption();

    QCommandLineOption inputOption("fileIn", "The input file <in>.", "in", QCoreApplication::applicationDirPath() + "/MNE-sample-data/MEG/sample/sample_audvis_raw.fif");
    QCommandLineOption outputOption("fileOut", "The output file <out>.", "out", QCoreApplication::applicationDirPath() + "/MNE-sample-data/MEG/sample/sample_audvis_filt_raw.fif");

    parser.addOption(inputOption);
    parser.addOption(outputOption);

    parser.process(a);

    // Init data loading
    QFile t_fileIn(parser.value(inputOption));
    QFile t_fileOut(parser.value(outputOption));

//    QFile t_fileIn("/Users/vorwerk/git/MNE-sample-data/EEG_data_001_raw_melon.fif");
//    QFile t_fileOut("/Users/vorwerk/git/MNE-sample-data/EEG_data_melon_bandpower_test.fif");

    FiffRawData raw(t_fileIn);

    // Set up the reading parameters to read the whole file at once
    fiff_int_t from = raw.first_samp;
    fiff_int_t to = raw.last_samp;

    // Read the data
    MatrixXd data;
    MatrixXd times;

    // Reading
    if(!raw.read_raw_segment(data, times, from, to)) {
        printf("error during read_raw_segment\n");
        return -1;
    }

    // Select channels and get indices

//    QStringList sUsedChannels;

/*    //FC3
    sUsedChannels.append("EEGFC3");
    sUsedChannels.append("EEGF3");
    sUsedChannels.append("EEGFC1");
    sUsedChannels.append("EEGC3");
    sUsedChannels.append("EEGFC5");

    //C3
    sUsedChannels.append("EEGC3");
    sUsedChannels.append("EEGFC3");
    sUsedChannels.append("EEGC1");
    sUsedChannels.append("EEGCP3");
    sUsedChannels.append("EEGC5");

    //CP3
    sUsedChannels.append("EEGCP3");
    sUsedChannels.append("EEGC3");
    sUsedChannels.append("EEGCP1");
    sUsedChannels.append("EEGP3");
    sUsedChannels.append("EEGCP5");

    //FC4
    sUsedChannels.append("EEGFC4");
    sUsedChannels.append("EEGF4");
    sUsedChannels.append("EEGFC2");
    sUsedChannels.append("EEGC4");
    sUsedChannels.append("EEGFC6");

    //C4
    sUsedChannels.append("EEGC4");
    sUsedChannels.append("EEGFC4");
    sUsedChannels.append("EEGC2");
    sUsedChannels.append("EEGCP4");
    sUsedChannels.append("EEGC6");

    //CP4
    sUsedChannels.append("EEGCP4");
    sUsedChannels.append("EEGC4");
    sUsedChannels.append("EEGCP2");
    sUsedChannels.append("EEGP4");
    sUsedChannels.append("EEGCP6");


    qDebug() << "sUsedChannels.length() " << sUsedChannels.length();

    VectorXi iUsedChannels(sUsedChannels.length());

   qDebug() << "iUsedChannels.size() " << iUsedChannels.size();

//    for(int i=0; i < raw.info.ch_names.length(); ++i)
//        qDebug() << raw.info.ch_names.at(i);

    for(int i=0; i < sUsedChannels.length(); ++i)
    {
        iUsedChannels[i] = raw.info.pick_channels(raw.info.ch_names,sUsedChannels.mid(i,1))[0];
        qDebug() << "Channel " << sUsedChannels.at(i) << " found at index " << iUsedChannels[i];
    }*/

    VectorXi iUsedChannels = raw.info.pick_types(false,true,false);
    for(int i=0; i < iUsedChannels.size(); ++i)
        qDebug() << raw.info.ch_names.at(iUsedChannels[i]);

    // Initialize filter settings
//    QString filter_name =  "Cosine_BPF";
//    FilterData::FilterType type = FilterData::BPF;
    double sFreq = raw.info.sfreq;
/*    double dCenterfreq = 10;
    double dBandwidth = 10;
    double dTransition = 1;

    RtFilter rtFilter;
    MatrixXd dataFiltered;*/

    // Set bandpower filter

    //int iBandPowerChannels = iUsedChannels.size()/5;
    int iBandPowerChannels = iUsedChannels.size();
    int iBandPowerBins = 11;

    double dUpdateIntervallLength = 0.05; // in [s]
    int iIntervallLengthFactor = 8; // multiples of update intervall

    int iUpdateIntervallLength = int(sFreq * dUpdateIntervallLength);
    int iIntervallLength = iIntervallLengthFactor * iUpdateIntervallLength;

    double dFreqOut = sFreq / double(iUpdateIntervallLength);

    qDebug() << "sFreq " << sFreq << " dFreqOut " << dFreqOut << " iUpdateInvervallLength " << iUpdateIntervallLength << " iIntervallLength " << iIntervallLength << "iBandPowerChannels " << iBandPowerChannels;

    int iOrderAR = 40;
    int iEvalsAR = 100;
    int iDetrendMethod = 2;
    double dFreqMin = 8;
    double dFreqMax = 30;


    // Only filter MEG and EEG channels
//    RowVectorXi picks_filter = raw.info.pick_types(true, true, false);

    // Get stim channels
    RowVectorXi stim_list = raw.info.pick_types(false, false, true);

    // Filtering
/*    printf("Filtering...");
    dataFiltered = rtFilter.filterData(data,
                                       type,
                                       dCenterfreq,
                                       dBandwidth,
                                       dTransition,
                                       sFreq,
                                       picks_filter);
    printf("[done]\n");*/

    // Get channels and rereference

    //MatrixXd datain = Eigen::MatrixXd::Zero(iUsedChannels.rows()/5,data.cols());
    MatrixXd datain = Eigen::MatrixXd::Zero(iUsedChannels.rows(),data.cols());

    qDebug() << "datain.rows() " << datain.rows() << " datain.cols() " << datain.cols();

//    for(int i=0; i<iUsedChannels.rows()/5; ++i)
//        datain.row(i) = data.row(iUsedChannels[i*5]) - 0.25 * (data.row(iUsedChannels[i*5]+1) + data.row(iUsedChannels[i*5]+2) + data.row(iUsedChannels[i*5]+3) + data.row(iUsedChannels[i*5]+4));

    for(int i=0; i<iUsedChannels.rows(); ++i)
        datain.row(i) = data.row(iUsedChannels[i]);
    RowVectorXd colavg = datain.colwise().sum();
    colavg /= iUsedChannels.size();
    for(int i=0; i<iUsedChannels.rows(); ++i)
        datain.row(i) -= colavg;

// Apply bandpower filter

    MatrixXd databp(iBandPowerChannels*iBandPowerBins,int(std::floor((datain.cols() - iIntervallLength)/iUpdateIntervallLength)) + 1);
    qDebug() << "databp.rows() " << databp.rows() << " databp.cols() " << databp.cols();

    qDebug() << "iDetrendMethod" << iDetrendMethod;

    for(int j=0; j<databp.cols(); j++)
    {
        MatrixXd datain_block = datain.block(0,j*iUpdateIntervallLength,datain.rows(),iIntervallLength);

        MatrixXd datain_noav = Spectral::detrendData(datain_block,iDetrendMethod);

        MatrixXcd ARSpectraWeights = Spectral::generateARSpectraWeights(dFreqMin/sFreq,dFreqMax/sFreq, iBandPowerBins, iEvalsAR, true);
        QVector<QPair<VectorXd, double>> ARCoeffs = Spectral::calculateARWeightsMEMMatrix(datain_noav,iOrderAR,true);
        QVector<VectorXd> matSpectrum = Spectral::psdFromARSpectra(ARCoeffs,ARSpectraWeights,sFreq,true);

        //stepwidth = (m_dFreqMax - m_dFreqMin)/(matSpectrum.at(0).size() - 1); //check where the evaluation points lie
        //for (int i=0; i < matSpectrum.length(); ++i)
        //        meanbandpower(0,0) += 1e10*bandpowerFromSpectrumEntries(matSpectrum.at(i),stepwidth)/matSpectrum.length();
        if(iBandPowerChannels < matSpectrum.length())
            qDebug() << "More channels selected than pre-defined! Only the first " << iBandPowerChannels << " are displayed!";
        for(int i=0; i<std::min(iBandPowerChannels,matSpectrum.length()); ++i){
            databp.block(i*iBandPowerBins,j,iBandPowerBins,1) = matSpectrum.at(i);
        }
        matSpectrum.clear();
    }

    // Create new fiffinfo

    FiffInfo out_fiff;// = raw.info;

    //out_fiff.filename = "t_fileOut";
    out_fiff.filename = "";

    out_fiff.nchan = iBandPowerChannels * iBandPowerBins;

    QList<FiffChInfo> fakeChList;
    FiffChInfo fakeCh;
    /*fakeCh.ch_name = "AR";
    fakeCh.kind = 502;
    fakeCh.range = -1;
    //fakeCh.setMinValue(0.5e-20);
    //fakeCh.setMaxValue(2e-20);
    fakeCh.unit = -1;
    fakeChList.append(fakeCh);
    fakeCh.ch_name = "FFT";
    fakeChList.append(fakeCh);*/
    fakeCh.ch_name = "BP ";
    fakeCh.kind = 502;
    fakeCh.range = -1;
    //fakeCh.setMinValue(0.5e-20);
    //fakeCh.setMaxValue(2e-20);
    fakeCh.unit = -1;

    QStringList fakeChNames;

    for(int i=0; i<iBandPowerChannels; ++i)
        for(int j=0; j<iBandPowerBins; ++j)
        {
            fakeCh.ch_name = QString("BP %1-%2").arg(raw.info.ch_names[iUsedChannels[i]]).arg(j);
            fakeChList.append(fakeCh);
            fakeChNames.append(fakeCh.ch_name);
        }

    for(int i=0; i<stim_list.size(); ++i)
    {
        fakeChList.append(raw.info.chs[stim_list[i]]);
        fakeChNames.append(raw.info.ch_names[stim_list[i]]);
    }

    out_fiff.chs = fakeChList;

    out_fiff.ch_names = fakeChNames;

    out_fiff.file_id = FIFFLIB::FiffId::new_file_id(); //check if necessary

    out_fiff.sfreq = dFreqOut;

    out_fiff.nchan = databp.rows() + stim_list.size();

    // Create output data

    MatrixXd dataout = MatrixXd::Zero(databp.rows() + stim_list.size(), databp.cols());
    int offset = databp.rows();

    for(int i=0; i<databp.rows(); ++i)
        dataout.row(i) = databp.row(i);

    for(int i=0; i<stim_list.size(); ++i)
    {
        RowVectorXd temprow = data.row(stim_list[i]);
        for(int j=0; j<databp.cols(); j++)
        {
            VectorXd tempsegment = temprow.segment(j*iUpdateIntervallLength,iUpdateIntervallLength);
            dataout((offset+i),j) = tempsegment.cwiseMax(0.0).maxCoeff();
        }
    }

    qDebug() << "dataout.rows()" << dataout.rows() << "dataout.cols()" << dataout.cols();

    // Init writing

    RowVectorXd cals;
    FiffStream::SPtr outfid = FiffStream::start_writing_raw(t_fileOut, out_fiff, cals);

    // Writing
    printf("Writing...");
    outfid->write_int(FIFF_FIRST_SAMPLE, &from);

    outfid->write_raw_buffer(dataout,cals);
    printf("[done]\n");

    outfid->finish_writing_raw();

    return 0;
}
