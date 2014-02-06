#ifdef DEBUG
  #define DEBUG_PRINT(...) do{ fprintf( stdout, __VA_ARGS__ ); } while( 0 )
  #define MSG_PRINT(...) do{ } while (0)
#else
  #define DEBUG_PRINT(...) do{ } while (0)
  #define MSG_PRINT(onoff, msg1, msg2) do { onoff > 0 ? fprintf(stderr, "\n%s\n\n%s\n\n", msg1, msg2) : fprintf(stderr, "\n%s\n\n", msg1); } while (0)
#endif
