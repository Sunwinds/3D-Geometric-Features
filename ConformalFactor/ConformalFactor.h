#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

class ConformalFactor: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "curvature.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Conformal Factor"; }
    QString description() { return "Conformal factor computation"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);
};

static QColor qtJetColorMap(double value, double min = 0.0, double max = 1.0){
    unsigned char rgb[3];
    unsigned char c1=144;
    float max4=(max-min)/4;
    value-=min;
    if(value==HUGE_VAL)
    {rgb[0]=rgb[1]=rgb[2]=255;}
    else if(value<0)
    {rgb[0]=rgb[1]=rgb[2]=0;}
    else if(value<max4)
    {rgb[0]=0;rgb[1]=0;rgb[2]=c1+(unsigned char)((255-c1)*value/max4);}
    else if(value<2*max4)
    {rgb[0]=0;rgb[1]=(unsigned char)(255*(value-max4)/max4);rgb[2]=255;}
    else if(value<3*max4)
    {rgb[0]=(unsigned char)(255*(value-2*max4)/max4);rgb[1]=255;rgb[2]=255-rgb[0];}
    else if(value<max)
    {rgb[0]=255;rgb[1]=(unsigned char)(255-255*(value-3*max4)/max4);rgb[2]=0;}
    else {rgb[0]=255;rgb[1]=rgb[2]=0;}
    return QColor(rgb[0],rgb[1],rgb[2]);
}
