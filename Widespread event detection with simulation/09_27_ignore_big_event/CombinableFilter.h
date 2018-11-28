#ifndef _COMBINABLEFILTER_H
#define _COMBINABLEFILTER_H

#pragma warning(disable:4996)

void CFS(long long* data, int numberOfAgents, char scheme, long long* monitorsArray, set<long long>* monitorEventSet, set<long long>** widespreadEventSet);
void CFPS(long long* data, int numberOfAgents, long long* monitorsArray, set<long long>* monitorEventSet, set<long long>** widespreadEventSet);

void CFS_Simulated(long long* data, int numberOfAgents, long long minEventSet, long long nCommonEvents, double Rint, double error, map<int, long long>* monitors, string scheme);
void CFSP_Simulated(long long* data, int numberOfAgents, long long minEventSet, long long nCommonEvents, double Rint, double error, map<int, long long>* monitors, string scheme);
void RawData_Simulated(long long* data, int numberOfAgents, long long averageEventSize, long long nCommonEvents, map<int, long long>* monitors);
void BDP_Simulated(long long* data, int numberOfAgents, long long averageEventSize, long long minEventSet, long long nCommonEvents, double Rint, double error, map<int, long long>* monitors, string scheme, int bitArraySize=0, int digest=0);
void SimulationforBitmapSize();
void SimulationCommunicationwrtRint();
void SimulationCommunicationwrtMonitors();
void SimulationAccuracyComparison();
void SimulationCommunicationwrtRintTheta();
void SimulationCommunicationwrtMonitorTheta();
void SimulationFalseEvents();
void SimulationFalseEventsBREDD();
void SimulationforFalsePositivewrtCommunicationOverhead();
#endif // !_PREPROCESS_H
