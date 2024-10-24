#include <QApplication>
#include <QMainWindow>
#include <QChartView>
#include <QScatterSeries>
#include <QValueAxis>
#include <QLineSeries>
#include <QFile>
#include <boost/json.hpp>
#include <vector>
#include <iostream>
#include "graphics.h"

using namespace QtCharts;
using namespace boost::json;

//Function to visulize the points and constraints
void visualizePoints(const std::vector<std::pair<double, double>>& points, const std::vector<std::pair<int, int>>& constraints) {
    //Creation of the main window for the visulization
    QMainWindow mainWindow;                                                      

    QScatterSeries* series = new QScatterSeries();  //Creation of a scatter series to print the points
    series->setMarkerSize(10);  //init the marker size to make points visible

    //Store the min and max of X and Y in order to set the limits of X and Y axis.
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();

    // Add points to the scatter series and calculate min/max values
    for (const auto& point : points) {
        double x = point.first;
        double y = point.second;

        series->append(x, y);

        // Find the min and max for both X,Y
        if (x < minX) minX = x;
        if (x > maxX) maxX = x;
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
    }

    //  Add a padding so the visulization to be cleaner (10% on both axis)
    double paddingX = (maxX - minX) * 0.1;
    double paddingY = (maxY - minY) * 0.1; 

    
    QLineSeries* lineSeries = new QLineSeries();// Creation of a line series for the constraints

    // Draw lines for each constraint
    for (const auto& constraint : constraints) {
        int point1_idx = constraint.first;
        int point2_idx = constraint.second;

        // Ensure the indices are valid
        if (point1_idx < points.size() && point2_idx < points.size()) {
            // Adding lines between the two points based on the constraint
            lineSeries->append(points[point1_idx].first, points[point1_idx].second);
            lineSeries->append(points[point2_idx].first, points[point2_idx].second);

            // Closing the line segment (re-adding the points)
            lineSeries->append(QPointF(points[point1_idx].first, points[point1_idx].second)); // Close the line segment
            lineSeries->append(QPointF(points[point2_idx].first, points[point2_idx].second)); // Close the line segment
        }
    }

    // Create a chart and add the series and constraints
    QChart* chart = new QChart();
    
    chart->addSeries(series); //For adding the points
    chart->addSeries(lineSeries); // For adding the constraints
    
    chart->setTitle("Point and Constraint Visualization");
    chart->createDefaultAxes();

    // Set the X axis
    QValueAxis* axisX = new QValueAxis();
    axisX->setTitleText("X Axis");
    axisX->setRange(minX - paddingX, maxX + paddingX); // Dynamically set X axis range with padding (based on the input points)

    // Set the Y axis
    QValueAxis* axisY = new QValueAxis();
    axisY->setTitleText("Y Axis");
    axisY->setRange(minY - paddingY, maxY + paddingY); // Dynamically set Y axis range with padding (based on the input points)

    //Attach the axes to the chart and associate them with the series
    chart->addAxis(axisX, Qt::AlignBottom);
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisX);
    series->attachAxis(axisY);
    lineSeries->attachAxis(axisX);
    lineSeries->attachAxis(axisY);

    // Creation of a chart view to the visulize the chart in the window
    QChartView* chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    mainWindow.setCentralWidget(chartView);
    mainWindow.resize(800, 600);
    mainWindow.show();

    // Start the loop for the application
    QApplication::exec();
}



