%extend casadi::DeserializerBase {
  %jscode %{
    DeserializerBase.prototype.unpack = function () {
      const name = SerializerBase.type_to_string(this.pop_type());
      return this["blind_unpack_" + name]();
    };
  %}
}
