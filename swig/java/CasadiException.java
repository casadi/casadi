package org.casadi;

/** Error raised by CasADi. */
public class CasadiException extends RuntimeException {
  private static final long serialVersionUID = 1L;

  public CasadiException(String message) {
    super(message);
  }

  public CasadiException(String message, Throwable cause) {
    super(message, cause);
  }
}
