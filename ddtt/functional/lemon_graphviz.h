#pragma once

#pragma warning(disable: 4267)
#define LEMON_ONLY_TEMPLATES
#include <lemon/list_graph.h>
#include <lemon/graph_to_eps.h>
#include <lemon/color.cc>
#include <lemon/bits/windows.cc>

#include <QStringList>

static QString digraphToGraphViz( lemon::ListDigraph & g, lemon::ListDigraph::NodeMap<QString> & nodeNames, QString caption = "Graph", QString subcaption = "" )
{
	QStringList out;
	out << "digraph G{\n";
	out << "\t" << "node [ fontcolor = black, color = white, style = filled ];" << "\n";

	// Nodes
	for (lemon::ListDigraph::NodeIt n(g); n != lemon::INVALID; ++n)
	{
		QString nodeName = nodeNames[n];
		QString shape = "rectangle";
		QString other = "";
		QString colorHex = "lightgrey";
		//colorHex.sprintf("#%02X%02X%02X", color.red(), color.green(), color.blue());

		out << "\t" << QString("%1 [label = \"%2\", color = \"%3\", shape = %4 %5];").arg(g.id(n)).arg( nodeName ).arg(colorHex).arg(shape).arg(other) << "\n";
	}

	// Edges
	for (lemon::ListDigraph::ArcIt a(g); a != lemon::INVALID; ++a)
	{
		QString color = "#000000";
		QString lable = "";
		QString lcolor = "#000000";

		out << "\t\"" << QString::number( g.id(g.source(a)) ) << "\" -> \"" << QString::number( g.id(g.target(a)) ) << "\"" 
			<< QString(" [color=\"%1\",label=\"%2\",fontcolor=\"%3\"] ").arg(color).arg(lable).arg(lcolor) << ";\n";
	}

	// Labels
	out << "label = \"\\n\\n" << caption << "\\n" << subcaption << "\"\n";
	out << "fontsize = 20;\n";

	out << "}\n";

	return out.join("");
}

void saveToGraphVizFile(QString str, QString filename = "graph")
{
	QFile file(filename + ".gv");
	if (!file.open(QFile::WriteOnly | QFile::Text)) return;
	QTextStream out(&file);
	out << str;
}
