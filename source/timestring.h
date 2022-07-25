#if !defined (_timestring_h_)
#define _timestring_h_
char *getTimeString();
void getISOdate(char *dateString, short maxLength);
void getISOtime(char *timeString, short maxLength);
void getISOdateTime(char *dateTimeString, short maxLength);
#endif
