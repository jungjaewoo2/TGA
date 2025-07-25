/*
 * coarseMCP3202.c
 * {azg] mods for Jaguar 
 * Thu Nov 30 21:38:27 PST 2017
 * run for 10 s only

rawMCP3202.c
Public Domain
2016-03-20

gcc -Wall -pthread -o coarseMCP3202  coarseMCP3202.c -lpigpio

This code shows how to bit bang SPI using DMA.

Using DMA to bit bang allows for two advantages

1) the time of the SPI transaction can be guaranteed to within a
   microsecond or so.

2) multiple devices of the same type can be read or written
  simultaneously.

This code shows how to read more than one MCP3202 at a time.

Each MCP3202 shares the SPI clock, MOSI, and slave select lines but has
a unique MISO line.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <pigpio.h>

#define SPI_SS 25 // GPIO for slave select.

#define ADCS 5    // Number of connected MCP3202.

#define BITS 12            // Bits per reading.
#define BX 6               // Bit position of data bit B11.
#define B0 (BX + BITS - 1) // Bit position of data bit B0.

#define MISO1 6   // ADC 1 MISO.
#define MISO2 26  //     2
#define MISO3 13  //     3
#define MISO4 23  //     4
#define MISO5 24  //     5

#define BUFFER 250       // Generally make this buffer as large as possible.

/* #define REPEAT_MICROS 40 // Reading every x microseconds.*/
#define REPEAT_MICROS 400 // Reading every x microseconds.

// #define SAMPLES 500000  // Number of samples to take,
// #define SAMPLES 100000  // Number of samples to take, ~ 40 sec 
// #define SAMPLES 62500  // Number of samples to take, ~ 25 sec 
#define SAMPLES 70000  // Number of samples to take, ~ 30 sec 

int MISO[ADCS]={MISO1, MISO2, MISO3, MISO4, MISO5};

rawSPI_t rawSPI =
{
   .clk     =  5, // GPIO for SPI clock.
   .mosi    = 12, // GPIO for SPI MOSI.
   .ss_pol  =  1, // Slave select resting level.
   .ss_us   =  1, // Wait 1 micro after asserting slave select.
   .clk_pol =  0, // Clock resting level.
   .clk_pha =  0, // 0 sample on first edge, 1 sample on second edge.
   .clk_us  =  1, // 2 clocks needed per bit so 500 kbps.
};

/*
   This function extracts the MISO bits for each ADC and
   collates them into a reading per ADC.
*/

void getReading(
   int adcs,  // Number of attached ADCs.
   int *MISO, // The GPIO connected to the ADCs data out.
   int OOL,   // Address of first OOL for this reading.
   int bytes, // Bytes between readings.
   int bits,  // Bits per reading.
   char *buf) 
{
   int i, a, p;
   uint32_t level;

   p = OOL;

   for (i=0; i<bits; i++)
   {
      level = rawWaveGetOut(p);

      for (a=0; a<adcs; a++)
      {
         putBitInBytes(i, buf+(bytes*a), level & (1<<MISO[a]));
      }

      p--;
   }
}


int main(int argc, char *argv[])
{
   int i, wid, offset;
   char buf[2];
   gpioPulse_t final[2];
   char rx[8];
   int sample;
   int val;
   int cb, botCB, topOOL, reading, now_reading;
   float cbs_per_reading;
   rawWaveInfo_t rwi;
   /* [azg] double start, end; */
   int pause;

   if (argc > 1) pause = atoi(argv[1]); else pause =0;

   if (gpioInitialise() < 0) return 1;

   // Need to set GPIO as outputs otherwise wave will have no effect.

   gpioSetMode(rawSPI.clk,  PI_OUTPUT);
   gpioSetMode(rawSPI.mosi, PI_OUTPUT);
   gpioSetMode(SPI_SS,      PI_OUTPUT);

   gpioWaveAddNew(); // Flush any old unused wave data.

   offset = 0;

   /*
   MCP3202 12-bit ADC 2 channels

   1  2  3  4  5  6   7   8  9  10 11 12 13 14 15 16 17
   SB SD OS MS NA B11 B10 B9 B8 B7 B6 B5 B4 B3 B2 B1 B0

   SB  1  1
   SD  1  0=differential 1=single
   OS  0  0=ch0, 1=ch1 (in single mode)
   MS  0  0=tx lsb first after tx msb first
   */

   /*
      Now construct lots of bit banged SPI reads.  Each ADC reading
      will be stored separately.  We need to ensure that the
      buffer is big enough to cope with any reasonable rescehdule.

      In practice make the buffer as big as you can.
   */

   for (i=0; i<BUFFER; i++)
   {
      buf[0] = 0xC0; // Start bit, single ended, channel 0.

      rawWaveAddSPI(&rawSPI, offset, SPI_SS, buf, 2, BX, B0, B0);

      /*
         REPEAT_MICROS must be more than the time taken to
         transmit the SPI message.
      */

      offset += REPEAT_MICROS;
   }

   // Force the same delay after the last reading.

   final[0].gpioOn = 0;
   final[0].gpioOff = 0;
   final[0].usDelay = offset;

   final[1].gpioOn = 0; // Need a dummy to force the final delay.
   final[1].gpioOff = 0;
   final[1].usDelay = 0;

   gpioWaveAddGeneric(2, final);

   wid = gpioWaveCreate(); // Create the wave from added data.

   if (wid < 0)
   {
      fprintf(stderr, "Can't create wave, %d too many?\n", BUFFER);
      return 1;
   }

   /*
      The wave resources are now assigned,  Get the number
      of control blocks (CBs) so we can calculate which reading
      is current when the program is running.
   */

   rwi = rawWaveInfo(wid);

   /* [azg]
    * printf("# cb %d-%d ool %d-%d del=%d ncb=%d nb=%d nt=%d\n",
      rwi.botCB, rwi.topCB, rwi.botOOL, rwi.topOOL, rwi.deleted,
      rwi.numCB,  rwi.numBOOL,  rwi.numTOOL);
      */

   printf("# matrix : d no_of_rows: %d no_of_cols: 1 REAL\n",SAMPLES);
   printf("# matrix columns 1 thru 5 \n");
   /*
      CBs are allocated from the bottom up.  As the wave is being
      transmitted the current CB will be between botCB and topCB
      inclusive.
   */

   botCB = rwi.botCB;

   /*
      Assume each reading uses the same number of CBs (which is
      true in this particular example).
   */

   cbs_per_reading = (float)rwi.numCB / (float)BUFFER;

   /*
    * [azg]
   printf("# cbs=%d per read=%.1f base=%d\n",
      rwi.numCB, cbs_per_reading, botCB);
      */

   /*
      OOL are allocated from the top down. There are BITS bits
      for each ADC reading and BUFFER ADC readings.  The readings
      will be stored in topOOL - 1 to topOOL - (BITS * BUFFER).
   */

   topOOL = rwi.topOOL;

   fprintf(stderr, "starting coarse test data collection...\n");

   if (pause) time_sleep(pause); // Give time to start a monitor.

   gpioWaveTxSend(wid, PI_WAVE_MODE_REPEAT);

   reading = 0;

   sample = 0;

   /* [azg] start = time_time(); */

   while (sample<SAMPLES)
   {
      // Which reading is current?

      cb = rawWaveCB() - botCB;

      now_reading = (float) cb / cbs_per_reading;

      // Loop gettting the fresh readings.

      while (now_reading != reading)
      {
         /*
            Each reading uses BITS OOL.  The position of this readings
            OOL are calculated relative to the waves top OOL.
         */
         getReading(
            ADCS, MISO, topOOL - ((reading%BUFFER)*BITS) - 1, 2, BITS, rx);

         /* [azg]  printf("%d", ++sample); */
	 ++sample; /* [azg] */

         /* for (i=0; i<ADCS; i++)*/
         for (i=4; i<ADCS; i++)
         {
            //   7   6  5  4  3  2  1  0 |  7  6  5  4  3  2  1  0
            // B11 B10 B9 B8 B7 B6 B5 B4 | B3 B2 B1 B0  X  X  X  X

            val = (rx[i*2]<<4) + (rx[(i*2)+1]>>4);
	    printf(" %d", val);
	     
         }

         printf("\n");

         if (++reading >= BUFFER) reading = 0;
      }
      usleep(1000);
   }

   /* [azg] end = time_time(); */

   /* [azg] printf("# %d samples in %.1f seconds (%.0f/s)\n",
      SAMPLES, end-start, (float)SAMPLES/(end-start)); */

   fprintf(stderr, "ending collection of coarse test data...\n");

   if (pause) time_sleep(pause);

   gpioTerminate();

   return 0;
}

