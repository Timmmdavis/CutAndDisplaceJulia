println("HelloREPL")

using Logging, IOLogging
logger = FileLogger(Dict(Logging.Info => "info.log", Logging.Error => "error.log"))
oldGlobal = global_logger(logger)
@info "Hello World!"
