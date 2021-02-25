package bhcd;

import bhcd.BasicDendrogram;
import bhcd.BasicTreeNode;
import bhcd.GraphAdjMatrix;
import bhcd.UniformDendrogram;
import blang.validation.ExactInvarianceTest;
import blang.validation.Instance;
import java.io.File;
import java.util.function.Function;
import org.eclipse.xtext.xbase.lib.Functions.Function1;
import org.eclipse.xtext.xbase.lib.IterableExtensions;
import org.eclipse.xtext.xbase.lib.ObjectExtensions;
import org.eclipse.xtext.xbase.lib.Procedures.Procedure1;
import org.junit.Test;

@SuppressWarnings("all")
public class TestUniformDendrogram {
  @Test
  public void invariance() {
    File _file = new File("data/terrorists/terrorist_pairs.csv");
    final GraphAdjMatrix graph = GraphAdjMatrix.parseToGraphAdjMatrix(_file);
    final UniformDendrogram model = new UniformDendrogram.Builder().setDend(BasicDendrogram.initDend(graph)).build();
    final Function<UniformDendrogram, Double> _function = (UniformDendrogram it) -> {
      double _xifexpression = (double) 0;
      final Function1<BasicTreeNode, String> _function_1 = (BasicTreeNode it_1) -> {
        return it_1.toString();
      };
      boolean _contains = IterableExtensions.<String>toList(IterableExtensions.<BasicTreeNode, String>map(it.getDend().getTree().successors(it.getDend().getRoot()), _function_1)).contains("1");
      if (_contains) {
        _xifexpression = 1.0;
      } else {
        _xifexpression = 0.0;
      }
      return Double.valueOf(_xifexpression);
    };
    final Instance<UniformDendrogram> instance = new Instance<UniformDendrogram>(model, _function);
    final Function<UniformDendrogram, Double> _function_1 = (UniformDendrogram it) -> {
      double _xifexpression = (double) 0;
      int _countLRleaves = it.getDend().countLRleaves(it.getDend().getRoot(), Boolean.valueOf(true));
      int _countLRleaves_1 = it.getDend().countLRleaves(it.getDend().getRoot(), Boolean.valueOf(false));
      boolean _equals = (_countLRleaves == _countLRleaves_1);
      if (_equals) {
        _xifexpression = 1.0;
      } else {
        _xifexpression = 0.0;
      }
      return Double.valueOf(_xifexpression);
    };
    final Instance<UniformDendrogram> instance2 = new Instance<UniformDendrogram>(model, _function_1);
    ExactInvarianceTest _exactInvarianceTest = new ExactInvarianceTest();
    final Procedure1<ExactInvarianceTest> _function_2 = (ExactInvarianceTest it) -> {
      it.add(instance);
      it.add(instance2);
    };
    final ExactInvarianceTest test = ObjectExtensions.<ExactInvarianceTest>operator_doubleArrow(_exactInvarianceTest, _function_2);
    test.check();
  }
}
