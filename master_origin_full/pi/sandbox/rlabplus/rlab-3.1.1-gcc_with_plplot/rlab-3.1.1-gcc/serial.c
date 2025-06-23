//
// FILE:
// serial.c
//
// HISTORY:
// Uses segments of code from project libserial.
// Copyright (c) 1995, 1997 Linas Vepstas
// Released under the  conditions  of  the  GNU General
// Public License.
//

#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <termios.h>
#include <unistd.h>

// ============================================================

static int serial_SetDTR ( Rfile * rf, int state )
{
  if ( !rf->fileds_i )
    goto _exit_serial_SetDTR;

  int set, rval=-1;
  ioctl (rf->fileds_i, TIOCMGET, &set);
  switch (state)
  {
    case -1:
      if (set & TIOCM_DTR)
        rval = 0x01;
      else
        rval = 0x00;
      break;

    case 1:
      set |= TIOCM_DTR;
      ioctl (rf->fileds_i, TIOCMSET, &set);
      rval = 1;
      break;

    case 0:
      set &= ~TIOCM_DTR;
      ioctl (rf->fileds_i, TIOCMSET, &set);
      rval = 0;
      break;

    default:
      break;
  }

_exit_serial_SetDTR:

  return rval;
}

static int serial_SetRTS ( Rfile * rf, int state )
{
  if ( !rf->fileds_i )
    goto _exit_serial_SetRTS;

  int set, rval=-1;
  ioctl (rf->fileds_i, TIOCMGET, &set);
  switch (state)
  {
    case -1:
      if (set & TIOCM_RTS)
        rval = 0x01;
      else
        rval = 0x00;
      break;

    case 1:
      set |= TIOCM_RTS;
      ioctl (rf->fileds_i, TIOCMSET, &set);
      rval = 1;
      break;

    case 0:
      set &= ~TIOCM_RTS;
      ioctl (rf->fileds_i, TIOCMSET, &set);
      rval = 0;
      break;

    default:
      break;
  }

_exit_serial_SetRTS:

  return rval;
}

static int serial_SetupCommPort (
    int fd, int ibrate, int idatb, int iparity, int istop , int iflow,
    int iraw, int hupcl
                                )
{
  int retval;

  if ( !fd )
    return 1;

  struct termios tset;

  // flush any unwritten, unread data
  tcflush(fd, TCIOFLUSH);
  sleep(1);

  // get current configuration of the port
  retval = tcgetattr(fd, &tset);
  if (-1 == retval)
    return errno;

  // Setup the terminal for ordinary I/O
  // e.g. as an ordinary modem or serial port printer.
  tset.c_cflag  = CLOCAL|CREAD|ibrate;
  tset.c_cflag &= ~CSIZE;

  //
  // set data format
  //
  // no. of data bits
  switch (idatb)
  {
    case 8:
      tset.c_cflag |= CS8;
      //tset.c_iflag &= ~ISTRIP;
      break;
    case 7:
      tset.c_cflag |= CS7;
      tset.c_iflag |= ISTRIP;
      break;
    case 6:
      tset.c_cflag |= CS6;
      tset.c_iflag |= ISTRIP;
      break;
    case 5:
      tset.c_cflag |= CS5;
      tset.c_iflag |= ISTRIP;
      break;
  }
  // parity
  switch (iparity)
  {
    case 0:
      // none
      tset.c_cflag &= ~PARENB;
      tset.c_iflag &= ~INPCK;
//       tset.c_iflag &=  IGNPAR; // was this wrong?
      tset.c_iflag |=  IGNPAR;
      break;
    case 1:
      // odd
      tset.c_cflag |=  PARENB;
      tset.c_cflag |=  PARODD;
      tset.c_iflag |=  INPCK;
      break;
    case 2:
      // even
      tset.c_cflag |=  PARENB;
      tset.c_cflag &= ~PARODD;
      tset.c_iflag |=  INPCK;
      break;
  }
  // no. of stop bits
  if (istop == 2)
    tset.c_cflag |=  CSTOPB;
  else
    tset.c_cflag &= ~CSTOPB;

  //
  // set flow control
  //
  tset.c_cflag     &= ~CRTSCTS;               // disable hardware
  tset.c_iflag     &= ~(IXON|IXOFF|IXANY);    // disable xon|xoff
  if (iflow == 2)
  {
    // enable hardware
    tset.c_cflag     |=  CRTSCTS;
    tset.c_cc[VSTART] = _POSIX_VDISABLE;
    tset.c_cc[VSTOP]  = _POSIX_VDISABLE;
  }
  else if (iflow == 1)
  {
    // enable xon|xoff
    tset.c_iflag     |=  IXON|IXOFF|IXANY;
  }

  if (hupcl)
    tset.c_cflag |=  HUPCL;
  else
    tset.c_cflag &= ~HUPCL;

  //
  // set defaults:
  //
  // no delay on carriage return, backspace, tab, etc.
  tset.c_oflag |= OPOST;
  tset.c_oflag |= NL0|CR0|TAB0|BS0|VT0|FF0;
  tset.c_oflag &= ~ (NLDLY|CRDLY|TABDLY|BSDLY|VTDLY|FFDLY);
  if (iraw)
  {
    // use raw input
    tset.c_lflag &= ~(ECHO | ECHONL | ICANON | ISIG | IEXTEN);

    // Default config: ignore break, do not ignore parity
    tset.c_iflag |= IGNBRK;

    // reset all output flags
    tset.c_oflag &= ~OPOST;

    // new:
//     tset.c_cflag &= ~(CSIZE | PARENB);
//     tset.c_cflag |= CS8;
    // new:
//     tset.c_iflag &= ~(IGNBRK | BRKINT | PARMRK | ISTRIP | INLCR | IGNCR | ICRNL | IXON);
  }

  //  default c_cc flags
  tset.c_cc[VEOF]     = CEOF;
  tset.c_cc[VQUIT]    = CQUIT;
  tset.c_cc[VKILL]    = CKILL;
  tset.c_cc[VEOL]     = CEOL;
  tset.c_cc[VSTART]   = CSTART;
  tset.c_cc[VSTOP]    = CSTOP;
  tset.c_cc[VSUSP]    = CSUSP;
  tset.c_cc[VREPRINT] = CRPRNT;
  tset.c_cc[VDISCARD] = CFLUSH;
  tset.c_cc[VWERASE]  = CWERASE;
  tset.c_cc[VLNEXT]   = CLNEXT;
  tset.c_cc[VMIN]     = 1;  // one character in the buffer
  tset.c_cc[VTIME]    = 1;  // do not really wait

  struct termios new_termios;

  retval = tcsetattr(fd, TCSANOW, &tset); // set new values
  if (-1 == retval)
    return errno;

  retval = 0;

  int l, k=0, M=500;
  unsigned char *t1= (unsigned char *)&tset, *t2=(unsigned char *)&new_termios;
  while (k < M)
  {
    k++;
    tcgetattr( fd, &new_termios );  // check new values
    if ( memcmp(t1, t2, sizeof(tset)) == 0 )
    {
      break;
    }
    usleep(1000);
    if (k == M)
    {
      fprintf(stderr, "Warning: Unable to set all serial port attributes!\n");
      retval = 2;
      if (RLABPLUS_SERIAL_DEBUG)
        for (l=0; l<sizeof(tset); l++)
        {
          if (t1[l]!=t2[l])
            fprintf(stderr, "%i: tset[%i]=%x, new_termios[%i]=%x\n", (int) sizeof(tset),
                    (int) l, (int) t1[l], (int) l, (int) t2[l]);
        }
    }
  }

  // flush: first attempt!
  tcsetattr(fd, TCSAFLUSH, &tset);

  sleep(2); //required to make flush work, for some reason
  tcflush(fd,TCIOFLUSH);

  return retval;
}

// ============================================================

static int int2baud ( int spd )
{
  int ibrate;

  // Always round up, so that we can talk to the
  // serial port at a speed equal or faster than
  // requested.
  // Flow control at excessive speeds is up to the user.
  if (460800 < spd) ibrate = 0;
  else if (230400 < spd) ibrate = B460800;
  else if (115200 < spd) ibrate = B230400;
  else if (57600 < spd)  ibrate = B115200;
  else if (38400 < spd)  ibrate = B57600;
  else if (19200 < spd)  ibrate = B38400;
  else if (14400 < spd)  ibrate = B19200;
  else if (4800 < spd)   ibrate = B9600;
  else if (2400 < spd)   ibrate = B4800;
  else if (1800 < spd)   ibrate = B2400;
  else if (1200 < spd)   ibrate = B1800;
  else if (600 < spd)    ibrate = B1200;
  else if (300 < spd)    ibrate = B600;
  else if (0 < spd)      ibrate = B300;
  else
    ibrate = 0;

  return ibrate;
}

// ============================================================

static int serial_Close ( Rfile * rf )
{
  struct termios tset;

  if (!rf->fileds_i) return 0;
  if (0x0 == rf->fileds_f) return 0;

  if ( rf->fileds_i )
  {
    int retval = tcgetattr (rf->fileds_i, &tset);
    if (-1 == retval)
    {
      fclose (rf->fileds_f);
      rf->fileds_f = 0x0;
      rf->fileds_i = 0;
      return errno;
    }

     // set terminal not to hangup on close
    tset.c_cflag &= ~HUPCL;
    tcsetattr (rf->fileds_i, TCSADRAIN, &tset);
  }

  fclose (rf->fileds_f);
  rf->fileds_f = 0x0;
  rf->fileds_i = 0;

  return 0;
}

// ============================================================

static int serial_OpenDevice (Rfile * rf,
                              int ispeed, int idatb, int iparity, int istop,
                              int iflow, int iraw, int hupcl
                             )
{
  int j;
  char * real_name=0;

  // ignore 'serial://' in front of the port name
  if (!strncmp(rf->name, "serial://", 9))
    real_name = &rf->name[9];
  else
    real_name = rf->name;

  //
  // if the device is already open do not do anything
  //
  if (rf->fileds_i > -1)
  {
    fflush (rf->fileds_f);
    return 0;
  }

  // check to see if this file exists already.
  // If it does not, then the user is trying to
  // create an ordinary file.
  struct stat buf;
  int retval = stat (real_name, &buf);
  if (retval)
  {
    printf ("Attempting to create ordinary file %s\n", real_name);
    rf->fileds_i = -1;
    return 1;
  }
  else
  {
    // open the serial port, make sure that its not the controlling tty
    // serialport = fopen ("/dev/modem", "r+");
    rf->fileds_i = 0;
    rf->fileds_f = 0;
    rf->fileds_i = open (real_name, O_RDWR | O_NOCTTY | O_NDELAY);
    rf->fileds_f = fdopen (rf->fileds_i, "a+");

    // if the device is e.g. a modem, and the modem is being used,
    // then busy-wait until it is free.
    j = 0;
    while ((!rf->fileds_f) && (EBUSY == errno) && (j < 10))
    {
      printf ("Serial Port %s appears to be busy, will retry "
          "%i-th time in 1 second\n", real_name, ++j);
      sleep (1);
      rf->fileds_i = open   (real_name, O_RDWR | O_NOCTTY | O_NDELAY);
      rf->fileds_f = fdopen (rf->fileds_i, "a+");
    }
    if (!rf->fileds_f)
      return errno; // not able to initialize serial port
  }

  // do the following only if not an ordinary file.
  if (rf->fileds_i > -1)
  {
    // flush any garbage remaining on the port from previous operations.
    sleep(2); //required to make flush work, for some reason
    retval = tcflush (rf->fileds_i, TCIOFLUSH);
    if (-1 == retval)
    {
      if (ENOTTY == errno) rf->fileds_i = 0;
      return errno;
    }

    // configure the port
    int ibrate = int2baud(ispeed);  // set baudrate
    serial_SetupCommPort (rf->fileds_i, ibrate, idatb, iparity, istop , iflow, iraw, hupcl );

    // allow commmunications to block
    int flags = fcntl(rf->fileds_i, F_GETFL, 0);
    retval    = fcntl(rf->fileds_i, F_SETFL, flags & ~O_NDELAY);
    if (retval == -1)
      return errno;
  }

  // everything is OK
  return 0;
}

// ====================== END OF FILE ============================
